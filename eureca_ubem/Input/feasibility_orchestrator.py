"""
feasibility_orchestrator.py

Scenario-swarm feasibility analysis for PUBEM.

This module is intentionally technical-only:
- no market levers
- no supplier tariff optimization
- no building NPV adoption
- no game steps

It samples controlled intervention diffusion points, generates random spatial
realisations, runs the existing PUBEM scenario engine, extracts grid/DHN
technical indicators, aggregates the swarm, and performs convergence analysis
on a continuous stress indicator.

Expected project imports
------------------------
The module assumes it is run inside the same environment where eureca_pubem is
available and where scenario_process.analyze_intervention can create a RETROFIT
scenario from an intervention dictionary.
"""

from __future__ import annotations

import copy
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd

from eureca_pubem import scenario_process as sc
from eureca_pubem import dhn_costs as dc
from eureca_pubem import grid_costs as gc


ENV_ORDER = ["none", "shallow", "medium", "deep"]
HP_ORDER = ["hp_le", "hp_me", "hp_he"]


# -----------------------------------------------------------------------------
# Data containers
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class PathConfig:
    baseline_geojson: str
    weatherfile_path: str
    output_folder: str = "./feasibility_outputs"

    dhn_pipe_cost_json: Optional[str] = None
    grid_cable_cost_json: Optional[str] = None


@dataclass(frozen=True)
class InterventionPolicy:
    """
    Defines how a diffusion point is translated into building configurations.

    conflict_mode:
        - "disjoint_hp_dhn": HP and DHN adopters are forced to be disjoint.
        - "split_end_uses": if a building is selected by both HP and DHN,
          space heating is assigned to HP and DHW is assigned to DHN.
    """

    target_envelope: str = "deep"
    target_hp_source: str = "hp_he"
    target_dhn_source: str = "dhn"
    target_pv_type: str = "A"
    target_pv_percentage: float = 100.0

    apply_hp_to: str = "both"       # "sh", "dhw", or "both"
    apply_dhn_to: str = "both"      # "sh", "dhw", or "both"
    conflict_mode: str = "disjoint_hp_dhn"

    allow_envelope_downgrade: bool = False
    allow_pv_decrease: bool = False

    boiler_fuel_fallback: str = "gas"


@dataclass(frozen=True)
class TechnicalLimits:
    min_voltage_pu: float = 0.95
    max_voltage_pu: float = 1.05
    max_supply_apparent_power: float = 1_000_000.0

    max_supply_heat: float = 500_000.0
    max_pipe_pressure: float = 2_000_000.0
    max_dp_per_m: float = 200.0
    max_pump_power: float = 1_000.0
    required_supply_temperature: float = 60.0


@dataclass(frozen=True)
class CostConfig:
    dhn_area_type: str = "urban"
    grid_area_type: str = "town"

    dhn_assumptions: Mapping[str, float] = field(default_factory=lambda: {
        "replacement_factor_ground": 0.5,
        "removal_factor": 0.3,
    })
    grid_assumptions: Mapping[str, float] = field(default_factory=lambda: {
        "replacement_factor_ground": 1.0,
        "removal_factor": 0.2,
        "ground_share": 0.55,
        "rest_share": 0.45,
    })
    grid_cable_key: str = "name"


@dataclass(frozen=True)
class DiffusionPoint:
    pv: float = 0.0
    hp: float = 0.0
    dhn: float = 0.0
    envelope: float = 0.0

    def as_key(self) -> str:
        return (
            f"pv={self.pv:.3f}|hp={self.hp:.3f}|"
            f"dhn={self.dhn:.3f}|env={self.envelope:.3f}"
        )
    
@dataclass
class FeasibilityPaths:
    baseline_geojson: str
    weatherfile_path: str
    dhn_pipe_cost_json: str
    grid_cable_cost_json: str

@dataclass
class RealisationResult:
    diffusion_key: str
    realisation_index: int
    seed: int
    n_buildings: int

    feasible: bool
    grid_feasible: bool
    dhn_feasible: bool
    combined_stress: float
    grid_stress: float
    dhn_stress: float

    indicators: Dict[str, float]
    intervention_counts: Dict[str, int]
    grid_faults: Dict[str, int]
    dhn_faults: Dict[str, int]

    grid_reinforcement_cost: Optional[float] = None
    dhn_reinforcement_cost: Optional[float] = None
    total_reinforcement_cost: Optional[float] = None


@dataclass
class SwarmAggregate:
    diffusion_key: str
    diffusion: DiffusionPoint
    n_realisations: int
    feasibility_probability: float
    grid_feasibility_probability: float
    dhn_feasibility_probability: float
    summary: Dict[str, float]


@dataclass
class ConvergenceRecord:
    diffusion_key: str
    r: int
    indicator_name: str
    indicator_value: float
    relative_change: Optional[float]


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------


def _ensure_output_folder(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _to_1d(x: Any) -> np.ndarray:
    if x is None:
        return np.array([], dtype=float)
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        return a.reshape(1)
    return a.reshape(-1)


def _safe_nanmax(x: Any, default: float = 0.0) -> float:
    a = _to_1d(x)
    if a.size == 0 or np.all(~np.isfinite(a)):
        return default
    return float(np.nanmax(a))


def _safe_nanmin(x: Any, default: float = 0.0) -> float:
    a = _to_1d(x)
    if a.size == 0 or np.all(~np.isfinite(a)):
        return default
    return float(np.nanmin(a))


def _safe_sum(x: Any, default: float = 0.0) -> float:
    a = _to_1d(x)
    if a.size == 0 or np.all(~np.isfinite(a)):
        return default
    return float(np.nansum(a))


def _ratio(value: float, limit: float, default: float = 0.0) -> float:
    if limit is None or not np.isfinite(limit) or limit <= 0:
        return default
    if value is None or not np.isfinite(value):
        return default
    return float(value) / float(limit)


def _normalize_source(value: Any) -> str:
    if value is None:
        return ""
    s = str(value).strip().lower()
    if "dhn" in s:
        return "dhn"
    if "hp_le" in s:
        return "hp_le"
    if "hp_me" in s:
        return "hp_me"
    if "hp_he" in s:
        return "hp_he"
    if s.startswith("hp") or "split" in s:
        return "hp_he" if s == "hp" else s
    if "boiler" in s:
        return "boiler"
    return s


def _env_rank(value: Any) -> int:
    s = str(value).strip().lower()
    return ENV_ORDER.index(s) if s in ENV_ORDER else 0


def _is_hp(value: Any) -> bool:
    return _normalize_source(value).startswith("hp")


def _target_count(share: float, n: int) -> int:
    share = float(np.clip(share, 0.0, 1.0))
    return int(round(share * n))


def _choose(rng: np.random.Generator, ids: Sequence[Any], n: int) -> List[Any]:
    ids = list(ids)
    if n <= 0 or not ids:
        return []
    n = min(int(n), len(ids))
    return list(rng.choice(np.asarray(ids, dtype=object), size=n, replace=False))


# -----------------------------------------------------------------------------
# Diffusion grid and intervention allocation
# -----------------------------------------------------------------------------


def build_diffusion_grid(
    pv_levels: Sequence[float],
    hp_levels: Sequence[float],
    dhn_levels: Sequence[float],
    envelope_levels: Sequence[float],
) -> List[DiffusionPoint]:
    points: List[DiffusionPoint] = []
    for pv in pv_levels:
        for hp in hp_levels:
            for dhn in dhn_levels:
                for env in envelope_levels:
                    points.append(
                        DiffusionPoint(
                            pv=float(pv),
                            hp=float(hp),
                            dhn=float(dhn),
                            envelope=float(env),
                        )
                    )
    return points


def load_baseline_gdf(path: str | Path) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path)
    required = ["id", "EEdepth", "SHSource", "DHWsource", "PVType", "PVpercentage"]
    missing = [c for c in required if c not in gdf.columns]
    if missing:
        raise ValueError(f"baseline_geojson is missing required columns: {missing}")
    return gdf


def eligible_ids_for_envelope(gdf: gpd.GeoDataFrame, target_env: str, allow_downgrade: bool) -> List[Any]:
    if allow_downgrade:
        return gdf["id"].tolist()
    target_rank = _env_rank(target_env)
    out = []
    for _, row in gdf.iterrows():
        if _env_rank(row["EEdepth"]) < target_rank:
            out.append(row["id"])
    return out


def eligible_ids_for_pv(gdf: gpd.GeoDataFrame, target_percentage: float, allow_decrease: bool) -> List[Any]:
    if allow_decrease:
        return gdf["id"].tolist()
    current = pd.to_numeric(gdf["PVpercentage"], errors="coerce").fillna(0.0)
    return gdf.loc[current < float(target_percentage), "id"].tolist()


def generate_intervention_realisation(
    baseline_gdf: gpd.GeoDataFrame,
    diffusion: DiffusionPoint,
    policy: InterventionPolicy,
    seed: int,
) -> Tuple[Dict[Any, Dict[str, Tuple[Any, Any]]], Dict[str, int]]:
    """
    Converts one diffusion point into one random intervention dictionary.

    The returned dictionary matches scenario_process.analyze_intervention(), where
    each changed field is represented as (old_value, new_value).
    """
    rng = np.random.default_rng(seed)
    gdf = baseline_gdf.copy()
    building_ids = gdf["id"].tolist()
    n = len(building_ids)

    n_pv = _target_count(diffusion.pv, n)
    n_hp = _target_count(diffusion.hp, n)
    n_dhn = _target_count(diffusion.dhn, n)
    n_env = _target_count(diffusion.envelope, n)

    env_ids = set(_choose(
        rng,
        eligible_ids_for_envelope(gdf, policy.target_envelope, policy.allow_envelope_downgrade),
        n_env,
    ))
    pv_ids = set(_choose(
        rng,
        eligible_ids_for_pv(gdf, policy.target_pv_percentage, policy.allow_pv_decrease),
        n_pv,
    ))

    if policy.conflict_mode == "disjoint_hp_dhn":
        hp_ids = set(_choose(rng, building_ids, n_hp))
        dhn_pool = [bid for bid in building_ids if bid not in hp_ids]
        if n_dhn > len(dhn_pool):
            raise ValueError(
                "Cannot allocate disjoint HP and DHN adopters: "
                f"hp={n_hp}, dhn={n_dhn}, buildings={n}. "
                "Use conflict_mode='split_end_uses' or reduce shares."
            )
        dhn_ids = set(_choose(rng, dhn_pool, n_dhn))
    elif policy.conflict_mode == "split_end_uses":
        hp_ids = set(_choose(rng, building_ids, n_hp))
        dhn_ids = set(_choose(rng, building_ids, n_dhn))
    else:
        raise ValueError("conflict_mode must be 'disjoint_hp_dhn' or 'split_end_uses'")

    interventions: Dict[Any, Dict[str, Tuple[Any, Any]]] = {bid: {} for bid in building_ids}

    row_by_id = gdf.set_index("id")

    for bid in env_ids:
        old = row_by_id.loc[bid, "EEdepth"]
        new = policy.target_envelope
        if policy.allow_envelope_downgrade or _env_rank(new) >= _env_rank(old):
            if old != new:
                interventions[bid]["envelope"] = (old, new)

    for bid in pv_ids:
        old_type = row_by_id.loc[bid, "PVType"]
        old_perc = float(row_by_id.loc[bid, "PVpercentage"] or 0.0)
        if old_type != policy.target_pv_type:
            interventions[bid]["PVtype"] = (old_type, policy.target_pv_type)
        new_perc = float(policy.target_pv_percentage)
        if policy.allow_pv_decrease or new_perc >= old_perc:
            if abs(old_perc - new_perc) > 1e-9:
                interventions[bid]["PVpercentage"] = (old_perc, new_perc)

    def apply_source(bid: Any, field: str, target: str) -> None:
        old = row_by_id.loc[bid, field]
        if old != target:
            key = "sh_source" if field == "SHSource" else "dhw_source"
            interventions[bid][key] = (old, target)

    for bid in hp_ids:
        if policy.conflict_mode == "split_end_uses" and bid in dhn_ids:
            # Split case: HP for SH, DHN for DHW by default.
            apply_source(bid, "SHSource", policy.target_hp_source)
            continue
        if policy.apply_hp_to in ("sh", "both"):
            apply_source(bid, "SHSource", policy.target_hp_source)
        if policy.apply_hp_to in ("dhw", "both"):
            apply_source(bid, "DHWsource", policy.target_hp_source)

    for bid in dhn_ids:
        if policy.conflict_mode == "split_end_uses" and bid in hp_ids:
            apply_source(bid, "DHWsource", policy.target_dhn_source)
            continue
        if policy.apply_dhn_to in ("sh", "both"):
            apply_source(bid, "SHSource", policy.target_dhn_source)
        if policy.apply_dhn_to in ("dhw", "both"):
            apply_source(bid, "DHWsource", policy.target_dhn_source)

    interventions = {bid: changes for bid, changes in interventions.items() if changes}

    counts = {
        "pv": len(pv_ids),
        "hp": len(hp_ids),
        "dhn": len(dhn_ids),
        "envelope": len(env_ids),
        "changed_buildings": len(interventions),
    }
    return interventions, counts


# -----------------------------------------------------------------------------
# Scenario execution
# -----------------------------------------------------------------------------


def run_single_realisation(
    *,
    diffusion: DiffusionPoint,
    realisation_index: int,
    seed: int,
    baseline_gdf: gpd.GeoDataFrame,
    city_demand_states: Mapping[str, Mapping[str, Any]],
    baseline_scenario: Any,
    paths: PathConfig,
    policy: InterventionPolicy,
    limits: TechnicalLimits,
    cost_config: Optional[CostConfig] = None,
    keep_scenario: bool = False,
) -> Tuple[RealisationResult, Optional[Any]]:
    interventions, counts = generate_intervention_realisation(
        baseline_gdf=baseline_gdf,
        diffusion=diffusion,
        policy=policy,
        seed=seed,
    )

    scenario, building_info, intervention_dictionary = sc.analyze_intervention(
        baseline_geojson=paths.baseline_geojson,
        city=city_demand_states,
        baseline_scenario=baseline_scenario,
        weatherfile_path=paths.weatherfile_path,
        mode="dictionary",
        intervention_dictionary=interventions,
    )

    result = extract_realisation_result(
        scenario=scenario,
        diffusion=diffusion,
        realisation_index=realisation_index,
        seed=seed,
        n_buildings=len(baseline_gdf),
        intervention_counts=counts,
        limits=limits,
        paths=paths,
        cost_config=cost_config,
    )

    return result, scenario if keep_scenario else None


# -----------------------------------------------------------------------------
# Indicator extraction
# -----------------------------------------------------------------------------


def count_grid_faults(scenario: Any) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for finding in getattr(scenario, "Electrical_Diagnoses", []) or []:
        code = str(getattr(finding, "code", "unknown"))
        out[code] = out.get(code, 0) + 1
    return out


def count_dhn_faults(scenario: Any) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for finding in getattr(scenario, "District_Heating_Diagnoses", []) or []:
        code = str(getattr(finding, "code", "unknown"))
        out[code] = out.get(code, 0) + 1
    return out


def extract_grid_indicators(scenario: Any, limits: TechnicalLimits) -> Tuple[Dict[str, float], float, bool]:
    max_vpu = -np.inf
    min_vpu = np.inf
    max_current_ratio = 0.0
    max_supply_s_ratio = 0.0
    max_backflow_w = 0.0
    all_converged = True

    n_grids = 0
    n_lines_overloaded = 0
    n_nodes_overvoltage = 0
    n_nodes_undervoltage = 0

    for grid in getattr(scenario, "Electrical_Network", []) or []:
        n_grids += 1
        res = getattr(grid, "results", None)
        if res is not None:
            Vm_pu = getattr(res, "Vm_pu", None)
            if Vm_pu is None:
                V = np.asarray(getattr(res, "V", []), dtype=float)
                slack = 400.0
                Vm_pu = V / slack if V.size else np.array([])
            Vm_pu = np.asarray(Vm_pu, dtype=float)
            if Vm_pu.size:
                max_vpu = max(max_vpu, float(np.nanmax(Vm_pu)))
                min_vpu = min(min_vpu, float(np.nanmin(Vm_pu)))
                n_nodes_overvoltage += int(np.any(Vm_pu > limits.max_voltage_pu, axis=0).sum()) if Vm_pu.ndim == 2 else int(np.any(Vm_pu > limits.max_voltage_pu))
                n_nodes_undervoltage += int(np.any(Vm_pu < limits.min_voltage_pu, axis=0).sum()) if Vm_pu.ndim == 2 else int(np.any(Vm_pu < limits.min_voltage_pu))

            S_supply = getattr(res, "S_supply", None)
            if S_supply is not None:
                max_supply_s_ratio = max(
                    max_supply_s_ratio,
                    _ratio(_safe_nanmax(S_supply), limits.max_supply_apparent_power),
                )

            P_supply = getattr(res, "P_supply", None)
            if P_supply is not None:
                pmin = _safe_nanmin(P_supply)
                if pmin < 0:
                    max_backflow_w = max(max_backflow_w, abs(pmin))

            converged = getattr(res, "converged", None)
            if converged is not None:
                all_converged = all_converged and bool(np.all(np.asarray(converged, dtype=bool)))

        for line in getattr(grid, "lines", []) or []:
            I = getattr(line, "I", None)
            if I is None and res is not None:
                line_index = getattr(res, "line_index", {}) or {}
                line_id = int(getattr(line, "line_id", -999999))
                if line_id in line_index:
                    I_line = np.asarray(getattr(res, "I_line", []), dtype=float)
                    if I_line.ndim == 2:
                        I = I_line[:, line_index[line_id]]
            max_current = float(getattr(line, "max_current_a", 0.0) or 0.0)
            r = _ratio(_safe_nanmax(I), max_current)
            max_current_ratio = max(max_current_ratio, r)
            if r > 1.0:
                n_lines_overloaded += 1

    if not np.isfinite(max_vpu):
        max_vpu = 0.0
    if not np.isfinite(min_vpu):
        min_vpu = 0.0

    voltage_high_stress = _ratio(max_vpu, limits.max_voltage_pu)
    voltage_low_stress = _ratio(limits.min_voltage_pu, min_vpu) if min_vpu > 0 else 0.0
    convergence_stress = 1.1 if not all_converged else 0.0

    grid_stress = max(
        voltage_high_stress,
        voltage_low_stress,
        max_current_ratio,
        max_supply_s_ratio,
        convergence_stress,
    )

    indicators = {
        "n_grids": float(n_grids),
        "grid_max_voltage_pu": float(max_vpu),
        "grid_min_voltage_pu": float(min_vpu),
        "grid_max_current_ratio": float(max_current_ratio),
        "grid_max_supply_s_ratio": float(max_supply_s_ratio),
        "grid_max_backflow_w": float(max_backflow_w),
        "grid_nodes_overvoltage": float(n_nodes_overvoltage),
        "grid_nodes_undervoltage": float(n_nodes_undervoltage),
        "grid_lines_overloaded": float(n_lines_overloaded),
        "grid_all_converged": float(all_converged),
        "grid_stress": float(grid_stress),
    }
    return indicators, float(grid_stress), bool(grid_stress <= 1.0)


def extract_dhn_indicators(scenario: Any, limits: TechnicalLimits) -> Tuple[Dict[str, float], float, bool]:
    n_systems = 0
    max_heat_ratio = 0.0
    max_pressure_ratio = 0.0
    max_dp_per_m_ratio = 0.0
    max_pump_ratio = 0.0
    supply_temperature_low_ratio = 0.0

    max_heat_injected = 0.0
    max_pressure_pa = 0.0
    max_dp_per_m = 0.0
    max_pump_power = 0.0
    min_supply_temp = np.inf
    total_heat_loss = 0.0

    for system in (getattr(scenario, "District_Heating_Systems", {}) or {}).values():
        n_systems += 1
        heat_loss = getattr(system, "hourly_heat_loss", None)
        total_heat_loss += _safe_sum(heat_loss)

        for node in getattr(system, "nodes", []) or []:
            heat_inj = getattr(node, "heat_injected", None)
            pump = getattr(node, "pump_power", None)
            sp = getattr(node, "supply_pressure", None)
            rp = getattr(node, "return_pressure", None)
            st = getattr(node, "supply_temperature", None)

            max_heat_injected = max(max_heat_injected, _safe_nanmax(heat_inj))
            max_pump_power = max(max_pump_power, _safe_nanmax(pump))
            max_pressure_pa = max(max_pressure_pa, _safe_nanmax(sp), _safe_nanmax(rp))
            if st is not None:
                min_supply_temp = min(min_supply_temp, _safe_nanmin(st, default=np.inf))

        for line in getattr(system, "lines", []) or []:
            L = float(getattr(line, "length", 0.0) or 0.0)
            if L <= 0:
                continue
            dp_s = _safe_nanmax(np.abs(_to_1d(getattr(line, "supply_pressure_drop", None))))
            dp_r = _safe_nanmax(np.abs(_to_1d(getattr(line, "return_pressure_drop", None))))
            max_dp_per_m = max(max_dp_per_m, dp_s / L, dp_r / L)
            max_line_p = float(getattr(line, "max_pressure", 0.0) or 0.0)
            if max_line_p > 0:
                max_pressure_ratio = max(max_pressure_ratio, _ratio(max_pressure_pa, max_line_p))

    if not np.isfinite(min_supply_temp):
        min_supply_temp = 0.0

    max_heat_ratio = _ratio(max_heat_injected, limits.max_supply_heat)
    max_pressure_ratio = max(max_pressure_ratio, _ratio(max_pressure_pa, limits.max_pipe_pressure))
    max_dp_per_m_ratio = _ratio(max_dp_per_m, limits.max_dp_per_m)
    max_pump_ratio = _ratio(max_pump_power, limits.max_pump_power)
    if min_supply_temp > 0:
        supply_temperature_low_ratio = _ratio(limits.required_supply_temperature, min_supply_temp)

    dhn_stress = max(
        max_heat_ratio,
        max_pressure_ratio,
        max_dp_per_m_ratio,
        max_pump_ratio,
        supply_temperature_low_ratio,
    )

    indicators = {
        "n_dhn_systems": float(n_systems),
        "dhn_max_heat_injected": float(max_heat_injected),
        "dhn_max_heat_ratio": float(max_heat_ratio),
        "dhn_max_pressure_pa": float(max_pressure_pa),
        "dhn_max_pressure_ratio": float(max_pressure_ratio),
        "dhn_max_dp_per_m": float(max_dp_per_m),
        "dhn_max_dp_per_m_ratio": float(max_dp_per_m_ratio),
        "dhn_max_pump_power": float(max_pump_power),
        "dhn_max_pump_ratio": float(max_pump_ratio),
        "dhn_min_supply_temperature": float(min_supply_temp),
        "dhn_supply_temperature_low_ratio": float(supply_temperature_low_ratio),
        "dhn_total_heat_loss": float(total_heat_loss),
        "dhn_stress": float(dhn_stress),
    }
    return indicators, float(dhn_stress), bool(dhn_stress <= 1.0)


def compute_reinforcement_costs(
    scenario: Any,
    paths: PathConfig,
    cost_config: Optional[CostConfig],
) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if cost_config is None:
        return None, None, None

    dhn_cost = None
    grid_cost = None

    if paths.dhn_pipe_cost_json is not None and hasattr(scenario, "dhn_pipe_changes"):
        pipe_table = dc.build_cost_table(paths.dhn_pipe_cost_json)
        dhn_cost = float(dc.compute_dhn_cost(
            dhn_pipe_changes=scenario.dhn_pipe_changes,
            area_type=cost_config.dhn_area_type,
            assumptions=dict(cost_config.dhn_assumptions),
            pipe_json=pipe_table,
        ))

    if paths.grid_cable_cost_json is not None and hasattr(scenario, "grid_line_changes"):
        grid_cost = float(gc.compute_grid_cost(
            grid_line_changes=scenario.grid_line_changes,
            area_type=cost_config.grid_area_type,
            assumptions=dict(cost_config.grid_assumptions),
            cable_json=paths.grid_cable_cost_json,
            cable_key=cost_config.grid_cable_key,
        ))

    if dhn_cost is None and grid_cost is None:
        total = None
    else:
        total = float((dhn_cost or 0.0) + (grid_cost or 0.0))
    return grid_cost, dhn_cost, total


def extract_realisation_result(
    *,
    scenario: Any,
    diffusion: DiffusionPoint,
    realisation_index: int,
    seed: int,
    n_buildings: int,
    intervention_counts: Dict[str, int],
    limits: TechnicalLimits,
    paths: PathConfig,
    cost_config: Optional[CostConfig] = None,
) -> RealisationResult:
    grid_ind, grid_stress, grid_feasible = extract_grid_indicators(scenario, limits)
    dhn_ind, dhn_stress, dhn_feasible = extract_dhn_indicators(scenario, limits)

    indicators = {**grid_ind, **dhn_ind}
    combined_stress = max(float(grid_stress), float(dhn_stress))
    feasible = bool(grid_feasible and dhn_feasible)

    grid_cost, dhn_cost, total_cost = compute_reinforcement_costs(scenario, paths, cost_config)

    return RealisationResult(
        diffusion_key=diffusion.as_key(),
        realisation_index=realisation_index,
        seed=seed,
        n_buildings=n_buildings,
        feasible=feasible,
        grid_feasible=bool(grid_feasible),
        dhn_feasible=bool(dhn_feasible),
        combined_stress=float(combined_stress),
        grid_stress=float(grid_stress),
        dhn_stress=float(dhn_stress),
        indicators=indicators,
        intervention_counts=intervention_counts,
        grid_faults=count_grid_faults(scenario),
        dhn_faults=count_dhn_faults(scenario),
        grid_reinforcement_cost=grid_cost,
        dhn_reinforcement_cost=dhn_cost,
        total_reinforcement_cost=total_cost,
    )


# -----------------------------------------------------------------------------
# Swarm execution and aggregation
# -----------------------------------------------------------------------------


def result_to_flat_dict(result: RealisationResult) -> Dict[str, Any]:
    row = {
        "diffusion_key": result.diffusion_key,
        "realisation_index": result.realisation_index,
        "seed": result.seed,
        "n_buildings": result.n_buildings,
        "feasible": result.feasible,
        "grid_feasible": result.grid_feasible,
        "dhn_feasible": result.dhn_feasible,
        "combined_stress": result.combined_stress,
        "grid_stress": result.grid_stress,
        "dhn_stress": result.dhn_stress,
        "grid_reinforcement_cost": result.grid_reinforcement_cost,
        "dhn_reinforcement_cost": result.dhn_reinforcement_cost,
        "total_reinforcement_cost": result.total_reinforcement_cost,
    }
    row.update({f"count_{k}": v for k, v in result.intervention_counts.items()})
    row.update(result.indicators)
    for k, v in result.grid_faults.items():
        row[f"grid_fault_{k}"] = v
    for k, v in result.dhn_faults.items():
        row[f"dhn_fault_{k}"] = v
    return row

def is_bool_like_series(series: pd.Series) -> bool:
    """
    Robust boolean detector.

    Avoids calling unique() on columns containing dict/list objects.
    """

    non_null = series.dropna()

    if non_null.empty:
        return False

    if pd.api.types.is_bool_dtype(non_null):
        return True

    # Only check object columns if all values are scalar bool-like.
    allowed = {True, False, 0, 1, "True", "False", "true", "false"}

    for value in non_null:
        # Skip dicts/lists/sets/tuples immediately.
        if isinstance(value, (dict, list, set, tuple)):
            return False

        if value not in allowed:
            return False

    return True


def aggregate_results(diffusion, results):
    """
    Aggregate realisation-level results for one diffusion point.

    Numeric columns:
        mean, std, p05, p50, p95, min, max

    Boolean columns:
        probability/share and count only

    This avoids pandas/numpy crashing when quantile is applied to booleans.
    """

    rows = []

    for result in results:
        if hasattr(result, "__dict__"):
            rows.append(result.__dict__)
        elif isinstance(result, dict):
            rows.append(result)
        else:
            raise TypeError(f"Unsupported result type: {type(result)}")

    df = pd.DataFrame(rows)

    summary = {
        "pv": diffusion.pv,
        "hp": diffusion.hp,
        "dhn": diffusion.dhn,
        "envelope": diffusion.envelope,
        "n_realisations": len(df),
    }

    for col in df.columns:
        if col in {"scenario", "scenario_obj", "interventions", "allocation", "realisation_index", "seed"}:
            continue

        series = df[col]

        if is_bool_like_series(series):
            b = series.dropna().astype(bool)

            if b.empty:
                continue

            summary[f"{col}_probability"] = float(b.mean())
            summary[f"{col}_count"] = int(b.sum())
            continue

        values = pd.to_numeric(series, errors="coerce").dropna()

        if values.empty:
            continue

        summary[f"{col}_mean"] = float(values.mean())
        summary[f"{col}_std"] = float(values.std(ddof=0))
        summary[f"{col}_p05"] = float(values.quantile(0.05))
        summary[f"{col}_p50"] = float(values.quantile(0.50))
        summary[f"{col}_p95"] = float(values.quantile(0.95))
        summary[f"{col}_min"] = float(values.min())
        summary[f"{col}_max"] = float(values.max())

    return summary

def run_swarm_for_diffusion(
    *,
    diffusion: DiffusionPoint,
    n_realisations: int,
    baseline_gdf: gpd.GeoDataFrame,
    city_demand_states: Mapping[str, Mapping[str, Any]],
    baseline_scenario: Any,
    paths: PathConfig,
    policy: InterventionPolicy,
    limits: TechnicalLimits,
    cost_config: Optional[CostConfig] = None,
    base_seed: int = 1000,
    keep_scenarios: bool = False,
    verbose: bool = True,
) -> Tuple[List[RealisationResult], SwarmAggregate, Optional[List[Any]]]:
    results: List[RealisationResult] = []
    scenarios: List[Any] = []

    for r in range(int(n_realisations)):
        seed = int(base_seed + 10_000 * abs(hash(diffusion.as_key())) % 1_000_000 + r)
        if verbose:
            print(f"{diffusion.as_key()} | realisation {r + 1}/{n_realisations} | seed={seed}")

        result, scenario = run_single_realisation(
            diffusion=diffusion,
            realisation_index=r,
            seed=seed,
            baseline_gdf=baseline_gdf,
            city_demand_states=city_demand_states,
            baseline_scenario=baseline_scenario,
            paths=paths,
            policy=policy,
            limits=limits,
            cost_config=cost_config,
            keep_scenario=keep_scenarios,
        )
        results.append(result)
        if keep_scenarios:
            scenarios.append(scenario)

    aggregate = aggregate_results(diffusion, results)
    return results, aggregate, scenarios if keep_scenarios else None


def run_feasibility_swarm(
    *,
    diffusion_points: Sequence[DiffusionPoint],
    n_realisations: int,
    baseline_geojson: str,
    city_demand_states: Mapping[str, Mapping[str, Any]],
    baseline_scenario: Any,
    weatherfile_path: str,
    policy: Optional[InterventionPolicy] = None,
    limits: Optional[TechnicalLimits] = None,
    cost_config: Optional[CostConfig] = None,
    dhn_pipe_cost_json: Optional[str] = None,
    grid_cable_cost_json: Optional[str] = None,
    output_folder: str = "./feasibility_outputs",
    base_seed: int = 1000,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[SwarmAggregate]]:
    """
    Runs the full feasibility swarm over a list of diffusion points.

    Parameters
    ----------
    city_demand_states:
        The dictionary returned by the envelope simulations, usually with keys
        such as "none", "shallow", "medium", "deep".
    baseline_scenario:
        The already created baseline/design scenario object used by
        scenario_process.analyze_intervention().
    """
    out_dir = _ensure_output_folder(output_folder)

    paths = PathConfig(
        baseline_geojson=baseline_geojson,
        weatherfile_path=weatherfile_path,
        output_folder=str(out_dir),
        dhn_pipe_cost_json=dhn_pipe_cost_json,
        grid_cable_cost_json=grid_cable_cost_json,
    )
    policy = policy or InterventionPolicy()
    limits = limits or TechnicalLimits()

    baseline_gdf = load_baseline_gdf(baseline_geojson)

    all_results: List[RealisationResult] = []
    aggregates: List[SwarmAggregate] = []

    for p_idx, diffusion in enumerate(diffusion_points):
        if verbose:
            print(f"\nDiffusion point {p_idx + 1}/{len(diffusion_points)}: {diffusion.as_key()}")

        results, aggregate, _ = run_swarm_for_diffusion(
            diffusion=diffusion,
            n_realisations=n_realisations,
            baseline_gdf=baseline_gdf,
            city_demand_states=city_demand_states,
            baseline_scenario=baseline_scenario,
            paths=paths,
            policy=policy,
            limits=limits,
            cost_config=cost_config,
            base_seed=base_seed,
            keep_scenarios=False,
            verbose=verbose,
        )
        all_results.extend(results)
        aggregates.append(aggregate)

    realisations_df = pd.DataFrame([result_to_flat_dict(r) for r in all_results])
    aggregate_rows = []
    
    for a in aggregates:
        if isinstance(a, dict):
            aggregate_rows.append(a)
        else:
            aggregate_rows.append(
                {
                    "diffusion_key": a.diffusion_key,
                    "pv": a.diffusion.pv,
                    "hp": a.diffusion.hp,
                    "dhn": a.diffusion.dhn,
                    "envelope": a.diffusion.envelope,
                    "n_realisations": a.n_realisations,
                    "feasibility_probability": a.feasibility_probability,
                    "grid_feasibility_probability": a.grid_feasibility_probability,
                    "dhn_feasibility_probability": a.dhn_feasibility_probability,
                    **a.summary,
                }
            )
    
    aggregates_df = pd.DataFrame(aggregate_rows)

    realisations_df.to_csv(out_dir / "realisations.csv", index=False)
    aggregates_df.to_csv(out_dir / "swarm_aggregates.csv", index=False)

    metadata = {
        "n_diffusion_points": len(diffusion_points),
        "n_realisations_per_point": int(n_realisations),
        "policy": asdict(policy),
        "limits": asdict(limits),
        "cost_config": None if cost_config is None else asdict(cost_config),
    }
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    return realisations_df, aggregates_df, aggregates


# -----------------------------------------------------------------------------
# Convergence analysis
# -----------------------------------------------------------------------------


def convergence_indicator_from_results(
    results: Sequence[RealisationResult],
    indicator_name: str = "combined_stress_p95",
) -> float:
    df = pd.DataFrame([result_to_flat_dict(r) for r in results])
    if indicator_name == "combined_stress_mean":
        return float(df["combined_stress"].mean())
    if indicator_name == "combined_stress_p95":
        return float(df["combined_stress"].quantile(0.95))
    if indicator_name == "grid_current_ratio_p95":
        return float(df["grid_max_current_ratio"].quantile(0.95))
    if indicator_name == "voltage_low_p05":
        return float(df["grid_min_voltage_pu"].quantile(0.05))
    if indicator_name == "voltage_high_p95":
        return float(df["grid_max_voltage_pu"].quantile(0.95))
    if indicator_name == "dhn_stress_p95":
        return float(df["dhn_stress"].quantile(0.95))
    if indicator_name == "reinforcement_cost_mean":
        col = "total_reinforcement_cost"
        if col not in df or df[col].isna().all():
            return float("nan")
        return float(df[col].mean())
    raise ValueError(f"Unknown convergence indicator: {indicator_name}")


def run_convergence_analysis(
    *,
    diffusion,
    r_values,
    baseline_geojson,
    city_demand_states,
    baseline_scenario,
    weatherfile_path,
    policy,
    limits,
    cost_config,
    dhn_pipe_cost_json,
    grid_cable_cost_json,
    indicator_name="combined_stress_p95",
    output_folder=None,
    base_seed=5000,
    verbose=True,
    relative_tolerance=0.01,
    consecutive_required=1,
):
    """
    Cumulative convergence analysis.

    Important:
    This does NOT rerun from 1 at every checkpoint.

    Example:
        r_values = [5, 10, 20, 40, 80]

    It runs:
        1..5
        then 6..10
        then 11..20
        then 21..40
        then 41..80

    It stops early if the selected convergence indicator stabilizes.

    Convergence criterion:
        relative_change <= relative_tolerance

    If consecutive_required > 1, the condition must hold for several checkpoints.
    """

    from pathlib import Path
    import numpy as np
    import pandas as pd
    import geopandas as gpd

    if output_folder is not None:
        out_dir = Path(output_folder)
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        out_dir = None

    r_values = sorted([int(r) for r in r_values])

    if not r_values:
        raise ValueError("r_values cannot be empty.")

    baseline_gdf = gpd.read_file(baseline_geojson)

    paths = FeasibilityPaths(
        baseline_geojson=baseline_geojson,
        weatherfile_path=weatherfile_path,
        dhn_pipe_cost_json=dhn_pipe_cost_json,
        grid_cable_cost_json=grid_cable_cost_json,
    )

    cumulative_results = []
    rows = []

    previous_value = None
    previous_r = 0
    stable_count = 0
    selected_r = None

    for target_r in r_values:
        if target_r <= previous_r:
            continue

        if verbose:
            print("\n" + "=" * 100)
            print(f"CONVERGENCE CHECKPOINT")
            print(f"Diffusion: {diffusion.as_key()}")
            print(f"Running additional realisations: {previous_r + 1} to {target_r}")
            print("=" * 100)

        for realisation_index in range(previous_r, target_r):
            seed = base_seed + realisation_index

            if verbose:
                print(
                    f"{diffusion.as_key()} | "
                    f"realisation {realisation_index + 1}/{target_r} | "
                    f"seed={seed}"
                )

            result, _ = run_single_realisation(
                diffusion=diffusion,
                realisation_index=realisation_index,
                seed=seed,
                baseline_gdf=baseline_gdf,
                city_demand_states=city_demand_states,
                baseline_scenario=baseline_scenario,
                paths=paths,
                policy=policy,
                limits=limits,
                cost_config=cost_config,
                keep_scenario=False,
            )

            cumulative_results.append(result)

        value = convergence_indicator_from_results(
            cumulative_results,
            indicator_name=indicator_name,
        )

        if previous_value is None:
            absolute_change = np.nan
            relative_change = np.nan
        else:
            absolute_change = abs(value - previous_value)

            if not np.isfinite(previous_value) or abs(previous_value) < 1e-12:
                relative_change = absolute_change
            else:
                relative_change = absolute_change / abs(previous_value)

        if np.isfinite(relative_change) and relative_change <= relative_tolerance:
            stable_count += 1
        else:
            stable_count = 0

        converged_here = stable_count >= consecutive_required

        row = {
            "diffusion_key": diffusion.as_key(),
            "pv": diffusion.pv,
            "hp": diffusion.hp,
            "dhn": diffusion.dhn,
            "envelope": diffusion.envelope,
            "r": target_r,
            "n_realisations": len(cumulative_results),
            "indicator_name": indicator_name,
            "indicator_value": value,
            "previous_indicator_value": previous_value,
            "absolute_change": absolute_change,
            "relative_change": relative_change,
            "relative_tolerance": relative_tolerance,
            "stable_count": stable_count,
            "consecutive_required": consecutive_required,
            "converged": converged_here,
        }

        rows.append(row)

        if verbose:
            print("\n" + "-" * 100)
            print(f"Convergence summary for {diffusion.as_key()}")
            print(f"R = {target_r}")
            print(f"{indicator_name} = {value}")
            print(f"previous = {previous_value}")
            print(f"absolute_change = {absolute_change}")
            print(f"relative_change = {relative_change}")
            print(f"stable_count = {stable_count}/{consecutive_required}")
            print(f"converged = {converged_here}")
            print("-" * 100)

        if out_dir is not None:
            pd.DataFrame(rows).to_csv(
                out_dir / "convergence_progress.csv",
                index=False,
            )

        previous_value = value
        previous_r = target_r

        if converged_here:
            selected_r = target_r

            if verbose:
                print("\n" + "!" * 100)
                print(f"CONVERGED for {diffusion.as_key()}")
                print(f"Selected R = {selected_r}")
                print(f"Stopping early.")
                print("!" * 100)

            break

    convergence_df = pd.DataFrame(rows)

    if selected_r is None:
        if "converged" in convergence_df.columns:
            valid = convergence_df[convergence_df["converged"] == True]
            if not valid.empty:
                selected_r = int(valid["r"].iloc[0])

    convergence_df["selected_r"] = selected_r

    if out_dir is not None:
        convergence_df.to_csv(out_dir / "convergence.csv", index=False)

    return convergence_df

# -----------------------------------------------------------------------------
# Example wiring
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    """
    Example only. Do not run this block before replacing the paths and before
    creating/loading `city_demand_states` and `baseline_scenario` in your own
    project session.

    Typical workflow in your notebook/script:

        from feasibility_orchestrator import *

        points = build_diffusion_grid(
            pv_levels=[0.0, 0.25, 0.50, 0.75, 1.0],
            hp_levels=[0.0, 0.25, 0.50, 0.75, 1.0],
            dhn_levels=[0.0, 0.25, 0.50, 0.75, 1.0],
            envelope_levels=[0.0, 0.50, 1.0],
        )

        realisations_df, aggregates_df, aggregates = run_feasibility_swarm(
            diffusion_points=points,
            n_realisations=40,
            baseline_geojson=PATHS["baseline_geojson"],
            city_demand_states=mycity,
            baseline_scenario=baseline,
            weatherfile_path=PATHS["weather_abs"],
            dhn_pipe_cost_json=PATHS["dhn_pipes"],
            grid_cable_cost_json=PATHS["cables"],
            cost_config=CostConfig(),
            output_folder="./feasibility_outputs",
        )

        convergence_df = run_convergence_analysis(
            diffusion=DiffusionPoint(pv=0.75, hp=0.75, dhn=0.0, envelope=1.0),
            r_values=[5, 10, 20, 40, 80],
            baseline_geojson=PATHS["baseline_geojson"],
            city_demand_states=mycity,
            baseline_scenario=baseline,
            weatherfile_path=PATHS["weather_abs"],
            indicator_name="combined_stress_p95",
        )
    """
    raise SystemExit("Import this module and call run_feasibility_swarm().")
