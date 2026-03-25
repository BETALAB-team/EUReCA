from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any,  List
import numpy as np


@dataclass
class DHEvidence:
    fault_type: str
    where: Any
    timesteps: List[int]
    metric_max: float
    limit: float


@dataclass
class DHFinding:
    code: str
    severity: str
    evidence: List[DHEvidence] = field(default_factory=list)


@dataclass
class DHDiagnosis:
    findings: List[DHFinding] = field(default_factory=list)


def _get(obj, name, default=None):
    if hasattr(obj, name):
        return getattr(obj, name)
    if isinstance(obj, dict) and name in obj:
        return obj[name]
    return default


def _as1d(x) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        return a.reshape(1)
    return a.reshape(-1)


def _severity_from_fraction(frac: float) -> str:
    if frac >= 0.05:
        return "critical"
    if frac >= 0.01:
        return "warning"
    return "info"


def _tsat_water_C_from_pressure_pa(p_pa: np.ndarray) -> np.ndarray:
    p_bar = np.clip(_as1d(p_pa) / 100000.0, 0.5, 25.0)
    P = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0], dtype=float)
    T = np.array([81.3, 100.0, 120.2, 133.5, 151.9, 179.9, 198.3, 212.4, 223.9], dtype=float)
    return np.interp(p_bar, P, T)


def _line_pressure_avg_pa(system, line, t: int, which: str) -> float:
    cache = _get(system, "pressure_cache", default=None)
    if cache is None:
        raise ValueError("system.pressure_cache missing")

    id_to_pos = cache.get("id_to_pos", None)
    if id_to_pos is None:
        node_ids = [int(n.node_id) for n in system.nodes]
        id_to_pos = {nid: i for i, nid in enumerate(node_ids)}

    s_id = int(getattr(line, "start_node"))
    e_id = int(getattr(line, "end_node"))

    if s_id not in id_to_pos or e_id not in id_to_pos:
        raise ValueError(f"Line {getattr(line,'line_id','?')} references unknown nodes")

    s = int(id_to_pos[s_id])
    e = int(id_to_pos[e_id])

    if which == "supply":
        Ps = float(_as1d(system.nodes[s].supply_pressure)[t])
        Pe = float(_as1d(system.nodes[e].supply_pressure)[t])
        return 0.5 * (Ps + Pe)

    Prs = float(_as1d(system.nodes[s].return_pressure)[t])
    Pre = float(_as1d(system.nodes[e].return_pressure)[t])
    return 0.5 * (Prs + Pre)


def diagnose_dhn_faults(
    system,
    *,
    max_supply_heat: float,
    boiling_margin: float,
    max_pipe_pressure: float,
    max_dp_per_m: float,
    max_pump_power: float,
) -> DHDiagnosis:
    findings: List[DHFinding] = []

    supply_nodes = [n for n in system.nodes if getattr(n, "node_type", None) == "supply"]
    if not supply_nodes:
        raise ValueError("No supply node in subsystem")
    supply = supply_nodes[0]

    heat_inj = _get(supply, "heat_injected", default=None)
    if heat_inj is None:
        raise ValueError("supply.heat_injected missing (run your metrics builder)")
    heat_inj = _as1d(heat_inj)
    T = int(heat_inj.shape[0])

    pump_power = _get(supply, "pump_power", default=None)
    if pump_power is None:
        raise ValueError("supply.pump_power missing (run your metrics builder)")
    pump_power = _as1d(pump_power)
    T = min(T, int(pump_power.shape[0]))

    def add_finding(code: str, evs: List[DHEvidence]):
        if not evs:
            return
        all_ts = sum((e.timesteps for e in evs), [])
        frac = len(set(all_ts)) / max(T, 1)
        findings.append(DHFinding(code=code, severity=_severity_from_fraction(frac), evidence=evs))

    ts = np.arange(T, dtype=int)

    bad = ts[heat_inj[:T] > float(max_supply_heat)]
    if bad.size:
        add_finding(
            "heat_overflow",
            [DHEvidence("heat_overflow", "supply", bad.tolist(), float(np.nanmax(heat_inj[:T])), float(max_supply_heat))],
        )

    bad = ts[pump_power[:T] > float(max_pump_power)]
    if bad.size:
        add_finding(
            "pump_power_exceeded",
            [DHEvidence("pump_power_exceeded", "supply", bad.tolist(), float(np.nanmax(pump_power[:T])), float(max_pump_power))],
        )

    ev_overT: List[DHEvidence] = []
    ev_overP: List[DHEvidence] = []
    ev_dpdx: List[DHEvidence] = []

    for line in system.lines:
        lid = int(getattr(line, "line_id", -1))
        L = float(getattr(line, "length", np.nan))
        if not np.isfinite(L) or L <= 0:
            continue

        Ts_s = _get(line, "supply_start_temperature", default=None)
        Ts_e = _get(line, "supply_end_temperature", default=None)
        Tr_s = _get(line, "return_start_temperature", default=None)
        Tr_e = _get(line, "return_end_temperature", default=None)

        if Ts_s is None or Ts_e is None or Tr_s is None or Tr_e is None:
            raise ValueError(f"Line {lid} missing temperature arrays (run thermal solve)")

        Ts_s = _as1d(Ts_s)[:T]
        Ts_e = _as1d(Ts_e)[:T]
        Tr_s = _as1d(Tr_s)[:T]
        Tr_e = _as1d(Tr_e)[:T]

        Tm_sup = 0.5 * (Ts_s + Ts_e)
        Tm_ret = 0.5 * (Tr_s + Tr_e)

        p_avg_sup = np.array([_line_pressure_avg_pa(system, line, t, "supply") for t in range(T)], dtype=float)
        p_avg_ret = np.array([_line_pressure_avg_pa(system, line, t, "return") for t in range(T)], dtype=float)

        Tsat_sup = _tsat_water_C_from_pressure_pa(p_avg_sup)
        Tsat_ret = _tsat_water_C_from_pressure_pa(p_avg_ret)

        Tmax_sup = (1.0 - float(boiling_margin)) * Tsat_sup
        Tmax_ret = (1.0 - float(boiling_margin)) * Tsat_ret

        bad_sup = ts[Tm_sup > Tmax_sup]
        bad_ret = ts[Tm_ret > Tmax_ret]
        badT = np.unique(np.concatenate([bad_sup, bad_ret])) if (bad_sup.size or bad_ret.size) else np.array([], dtype=int)

        if badT.size:
            metric_max = float(np.nanmax(np.concatenate([Tm_sup, Tm_ret])))
            limit_min = float(np.nanmin(np.concatenate([Tmax_sup, Tmax_ret])))
            ev_overT.append(DHEvidence("pipe_overtemperature", lid, badT.tolist(), metric_max, limit_min))

        lp = getattr(line, "max_pressure", None)
        limit_p = float(lp) if (lp is not None and np.isfinite(float(lp)) and float(lp) > 0) else float(max_pipe_pressure)

        badP = ts[(p_avg_sup > limit_p) | (p_avg_ret > limit_p)]
        if badP.size:
            metric_max = float(np.nanmax(np.concatenate([p_avg_sup, p_avg_ret])))
            ev_overP.append(DHEvidence("pipe_overpressure", lid, badP.tolist(), metric_max, float(limit_p)))

        dps = _get(line, "supply_pressure_drop", default=None)
        dpr = _get(line, "return_pressure_drop", default=None)
        if dps is None or dpr is None:
            raise ValueError(f"Line {lid} missing pressure drop arrays (run hydraulics)")

        dps = _as1d(dps)[:T]
        dpr = _as1d(dpr)[:T]

        dpdx = np.maximum(np.abs(dps), np.abs(dpr)) / L
        badD = ts[dpdx > float(max_dp_per_m)]
        if badD.size:
            ev_dpdx.append(
                DHEvidence(
                    "excessive_pressure_gradient",
                    lid,
                    badD.tolist(),
                    float(np.nanmax(dpdx)),
                    float(max_dp_per_m),
                )
            )

    add_finding("pipe_overtemperature", ev_overT)
    add_finding("pipe_overpressure", ev_overP)
    add_finding("excessive_pressure_gradient", ev_dpdx)

    return DHDiagnosis(findings=findings)
