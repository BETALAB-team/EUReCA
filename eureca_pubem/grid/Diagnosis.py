from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Tuple, Any, Optional
import numpy as np


def _to_1d(arr) -> np.ndarray:
    return np.asarray(arr, dtype=float).reshape(-1)


class Domain(str, Enum):
    ELECTRICAL = "electrical"


class Quantity(str, Enum):
    VOLTAGE = "voltage"
    CURRENT = "current"
    APPARENT_POWER = "apparent_power"
    ACTIVE_POWER = "active_power"


class FaultKind(str, Enum):
    THRESHOLD_HIGH = "threshold_high"
    THRESHOLD_LOW = "threshold_low"
    DIRECTION_REVERSE = "direction_reverse"


class Scope(str, Enum):
    NODE = "node"
    LINE = "line"
    GRID = "grid"


class Severity(str, Enum):
    INFO = "info"
    WARNING = "warning"
    CRITICAL = "critical"


@dataclass(frozen=True)
class FaultType:
    code: str
    domain: Domain
    quantity: Quantity
    kind: FaultKind
    scope: Scope


FAULT_TYPES: Dict[str, FaultType] = {
    "overvoltage": FaultType(
        code="overvoltage",
        domain=Domain.ELECTRICAL,
        quantity=Quantity.VOLTAGE,
        kind=FaultKind.THRESHOLD_HIGH,
        scope=Scope.NODE,
    ),
    "undervoltage": FaultType(
        code="undervoltage",
        domain=Domain.ELECTRICAL,
        quantity=Quantity.VOLTAGE,
        kind=FaultKind.THRESHOLD_LOW,
        scope=Scope.NODE,
    ),
    "overcurrent": FaultType(
        code="overcurrent",
        domain=Domain.ELECTRICAL,
        quantity=Quantity.CURRENT,
        kind=FaultKind.THRESHOLD_HIGH,
        scope=Scope.LINE,
    ),
    "overflow": FaultType(
        code="overflow",
        domain=Domain.ELECTRICAL,
        quantity=Quantity.APPARENT_POWER,
        kind=FaultKind.THRESHOLD_HIGH,
        scope=Scope.GRID,
    ),
    "backflow": FaultType(
        code="backflow",
        domain=Domain.ELECTRICAL,
        quantity=Quantity.ACTIVE_POWER,
        kind=FaultKind.DIRECTION_REVERSE,
        scope=Scope.GRID,
    ),
}


@dataclass
class FaultEvent:
    ftype: FaultType
    where: int  # node_index, line_index, or -1 for grid-level
    timesteps: List[int]
    metric_min: Optional[float] = None
    metric_max: Optional[float] = None
    limit: Optional[float] = None
    severity: Severity = Severity.WARNING
    tags: Dict[str, Any] = field(default_factory=dict)


@dataclass
class GridFaults:
    events_by_type: Dict[str, List[FaultEvent]] = field(default_factory=dict)
    summary: Dict[str, Any] = field(default_factory=dict)

    def add_event(self, ev: FaultEvent) -> None:
        self.events_by_type.setdefault(ev.ftype.code, []).append(ev)


@dataclass
class FaultReport:
    grids: Dict[int, GridFaults] = field(default_factory=dict)


def _severity_from_fraction(frac: float) -> Severity:
    # tweak these thresholds however you like
    if frac >= 0.20:
        return Severity.CRITICAL
    if frac >= 0.02:
        return Severity.WARNING
    return Severity.INFO


def generate_grid_fault_report_object(
    gs,
    min_voltage: float,
    max_voltage: float,
    max_supply_apparent_power: float,
) -> FaultReport:
    report = FaultReport()
    
    if type(gs) == list:
        whiler = gs
    else: 
        whiler = gs.grids

    for gi, grid in enumerate(whiler):
        gf = GridFaults()

        # We’ll compute some context features automatically
        T = None

        # Voltage faults (buildings only)
        V_max = -np.inf
        V_min = np.inf
        ov_nodes = 0
        uv_nodes = 0

        for ni, node in enumerate(grid.nodes):
            if getattr(node, "node_type", None) != "building":
                continue

            V = _to_1d(node.V)
            if V.size == 0:
                continue
            if T is None:
                T = int(V.size)

            V_max = max(V_max, float(np.nanmax(V)))
            V_min = min(V_min, float(np.nanmin(V)))

            ov = np.where(V > max_voltage)[0]
            if ov.size:
                ov_nodes += 1
                frac = float(ov.size) / float(V.size)
                gf.add_event(
                    FaultEvent(
                        ftype=FAULT_TYPES["overvoltage"],
                        where=ni,
                        timesteps=ov.astype(int).tolist(),
                        metric_min=float(np.nanmin(V[ov])),
                        metric_max=float(np.nanmax(V[ov])),
                        limit=float(max_voltage),
                        severity=_severity_from_fraction(frac),
                        tags={
                            "node_type": "building",
                            "violation_fraction": frac,
                        },
                    )
                )

            uv = np.where(V < min_voltage)[0]
            if uv.size:
                uv_nodes += 1
                frac = float(uv.size) / float(V.size)
                gf.add_event(
                    FaultEvent(
                        ftype=FAULT_TYPES["undervoltage"],
                        where=ni,
                        timesteps=uv.astype(int).tolist(),
                        metric_min=float(np.nanmin(V[uv])),
                        metric_max=float(np.nanmax(V[uv])),
                        limit=float(min_voltage),
                        severity=_severity_from_fraction(frac),
                        tags={
                            "node_type": "building",
                            "violation_fraction": frac,
                        },
                    )
                )

        # Overcurrent faults
        oc_lines = 0
        for li, line in enumerate(grid.lines):
            I = _to_1d(line.I)
            if I.size == 0:
                continue
            if T is None:
                T = int(I.size)

            limit = float(line.max_current_a)
            oc = np.where(I > limit)[0]
            if oc.size:
                oc_lines += 1
                frac = float(oc.size) / float(I.size)
                gf.add_event(
                    FaultEvent(
                        ftype=FAULT_TYPES["overcurrent"],
                        where=li,
                        timesteps=oc.astype(int).tolist(),
                        metric_min=float(np.nanmin(I[oc])),
                        metric_max=float(np.nanmax(I[oc])),
                        limit=limit,
                        severity=_severity_from_fraction(frac),
                        tags={
                            "violation_fraction": frac,
                        },
                    )
                )

        # Overflow (grid-level, apparent power)
        S = _to_1d(grid.results.S_supply)
        if S.size:
            if T is None:
                T = int(S.size)
            of = np.where(S > float(max_supply_apparent_power))[0]
            if of.size:
                frac = float(of.size) / float(S.size)
                gf.add_event(
                    FaultEvent(
                        ftype=FAULT_TYPES["overflow"],
                        where=-1,
                        timesteps=of.astype(int).tolist(),
                        metric_min=float(np.nanmin(S[of])),
                        metric_max=float(np.nanmax(S[of])),
                        limit=float(max_supply_apparent_power),
                        severity=_severity_from_fraction(frac),
                        tags={"violation_fraction": frac},
                    )
                )

        # Backflow (grid-level, active power)
        P = _to_1d(grid.results.P_supply)
        if P.size:
            if T is None:
                T = int(P.size)
            bf = np.where(P < 0.0)[0]
            if bf.size:
                frac = float(bf.size) / float(P.size)
                gf.add_event(
                    FaultEvent(
                        ftype=FAULT_TYPES["backflow"],
                        where=-1,
                        timesteps=bf.astype(int).tolist(),
                        metric_min=float(np.nanmin(P[bf])),
                        metric_max=float(np.nanmax(P[bf])),
                        limit=0.0,
                        severity=_severity_from_fraction(frac),
                        tags={"violation_fraction": frac},
                    )
                )

        # Grid summary (taxonomy-friendly “context”)
        gf.summary = {
            "T": T,
            "V_max": None if not np.isfinite(V_max) else float(V_max),
            "V_min": None if not np.isfinite(V_min) else float(V_min),
            "building_nodes_with_overvoltage": int(ov_nodes),
            "building_nodes_with_undervoltage": int(uv_nodes),
            "lines_with_overcurrent": int(oc_lines),
            "fault_types_present": sorted(gf.events_by_type.keys()),
        }

        report.grids[gi] = gf

    return report



from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Any, Optional, Tuple
import numpy as np


class DiagnosisCode(str, Enum):
    # Voltage-related
    WIDESPREAD_UNDERVOLTAGE = "widespread_undervoltage"
    LOCAL_UNDERVOLTAGE_WITH_OVERCURRENT = "local_undervoltage_with_overcurrent"
    WIDESPREAD_OVERVOLTAGE = "widespread_overvoltage"
    LOCAL_OVERVOLTAGE_WITH_BACKFLOW = "local_overvoltage_with_backflow"

    # Thermal/asset loading proxy
    LINE_OVERLOAD = "line_overload"
    SUPPLY_CAPACITY_EXCEEDED = "supply_capacity_exceeded"
    REVERSE_POWER_FLOW = "reverse_power_flow"

    # Mixed / ambiguous
    MIXED_VOLTAGE_ISSUES = "mixed_voltage_issues"
    UNKNOWN_PATTERN = "unknown_pattern"


@dataclass
class Evidence:
    fault_type: str                # e.g. "undervoltage"
    where: int                     # node idx, line idx, -1 grid
    timesteps: List[int]
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TroubleshootStep:
    code: str                      # stable machine key, e.g. "check_tap_changer"
    description: str               # short human text
    data_needed: List[str] = field(default_factory=list)  # what extra inputs would improve confidence


@dataclass
class DiagnosisFinding:
    code: DiagnosisCode
    title: str
    confidence: float              # 0..1
    severity: str                  # "info"|"warning"|"critical"
    evidence: List[Evidence] = field(default_factory=list)
    likely_causes: List[str] = field(default_factory=list)     # short labels
    next_steps: List[TroubleshootStep] = field(default_factory=list)
    notes: Dict[str, Any] = field(default_factory=dict)
    sentence: Optional[str] = None  # optional human sentence


@dataclass
class DiagnosisReport:
    grid_index: int
    findings: List[DiagnosisFinding] = field(default_factory=list)
    summary: Dict[str, Any] = field(default_factory=dict)


# -------- helpers (no magic, just useful primitives) --------

def _events(gf, fault_code: str):
    return gf.events_by_type.get(fault_code, []) if gf else []

def _unique_where(events) -> List[int]:
    return sorted({int(ev.where) for ev in events})

def _all_timesteps(events) -> np.ndarray:
    if not events:
        return np.array([], dtype=int)
    ts = np.concatenate([np.asarray(ev.timesteps, dtype=int) for ev in events if ev.timesteps])
    if ts.size == 0:
        return ts
    ts = np.unique(ts)
    return ts

def _intersect_timesteps(a: List[int], b: List[int]) -> List[int]:
    if not a or not b:
        return []
    A = set(a)
    return sorted(A.intersection(b))

def _coverage_fraction(events, T: Optional[int]) -> float:
    if T is None or T <= 0:
        return 0.0
    ts = _all_timesteps(events)
    return float(ts.size) / float(T)

def _severity_from_conf(conf: float) -> str:
    if conf >= 0.75:
        return "critical"
    if conf >= 0.4:
        return "warning"
    return "info"

def _mk_evidence_from_event(ev, fault_code: str) -> Evidence:
    return Evidence(
        fault_type=fault_code,
        where=int(ev.where),
        timesteps=list(ev.timesteps),
        details={
            "limit": ev.limit,
            "metric_min": ev.metric_min,
            "metric_max": ev.metric_max,
            "severity": getattr(ev, "severity", None),
            "tags": dict(ev.tags) if getattr(ev, "tags", None) else {},
        },
    )


# -------- core diagnosis engine --------

def diagnose_grid_faults(
    gf,
    grid_index: int,
    widespread_node_fraction: float = 0.25,   # “widespread” if >= 25% of building nodes affected
    widespread_time_fraction: float = 0.05,   # “sustained” if >= 5% of year
    coincidence_min_steps: int = 5,           # ignore tiny coincidences
) -> DiagnosisReport:
    T = gf.summary.get("T") if gf else None

    ov = _events(gf, "overvoltage")
    uv = _events(gf, "undervoltage")
    oc = _events(gf, "overcurrent")
    of = _events(gf, "overflow")
    bf = _events(gf, "backflow")

    n_bldg = None
    # we only know counts from summary in this v1
    ov_nodes = int(gf.summary.get("building_nodes_with_overvoltage", 0))
    uv_nodes = int(gf.summary.get("building_nodes_with_undervoltage", 0))
    oc_lines = int(gf.summary.get("lines_with_overcurrent", 0))

    findings: List[DiagnosisFinding] = []

    # --- Supply capacity exceeded (overflow) ---
    if of:
        frac = _coverage_fraction(of, T)
        conf = min(1.0, 0.6 + 2.0 * frac)  # more time over limit -> more confidence
        evs = [_mk_evidence_from_event(e, "overflow") for e in of]
        findings.append(
            DiagnosisFinding(
                code=DiagnosisCode.SUPPLY_CAPACITY_EXCEEDED,
                title="Supply apparent power limit exceeded",
                confidence=conf,
                severity=_severity_from_conf(conf),
                evidence=evs,
                likely_causes=[
                    "peak demand exceeds substation rating",
                    "unmodeled simultaneity / peak coincidence",
                    "DER export increases apparent power",
                ],
                next_steps=[
                    TroubleshootStep(
                        code="check_substation_rating_and_limits",
                        description="Verify transformer/inverter apparent power ratings and configured limits (S_max).",
                        data_needed=["transformer rating", "inverter ratings", "protection settings"],
                    ),
                    TroubleshootStep(
                        code="inspect_peak_hours",
                        description="Identify peak hours and which loads/generators dominate those hours.",
                        data_needed=["per-building P/Q time series", "PV/export profiles"],
                    ),
                ],
                sentence="Supply capacity is exceeded at some timesteps; investigate peak loading and substation S-limits.",
                notes={"overflow_time_fraction": frac},
            )
        )

    # --- Reverse power flow (backflow) ---
    if bf:
        frac = _coverage_fraction(bf, T)
        conf = min(1.0, 0.55 + 2.0 * frac)
        evs = [_mk_evidence_from_event(e, "backflow") for e in bf]
        findings.append(
            DiagnosisFinding(
                code=DiagnosisCode.REVERSE_POWER_FLOW,
                title="Reverse active power flow at supply",
                confidence=conf,
                severity=_severity_from_conf(conf),
                evidence=evs,
                likely_causes=[
                    "distributed generation export (PV) exceeds local demand",
                    "control policy allows export",
                ],
                next_steps=[
                    TroubleshootStep(
                        code="confirm_export_policy",
                        description="Check whether export is allowed/expected (tariffs, protection, control policy).",
                        data_needed=["export limits", "control policy", "tariff rules"],
                    ),
                    TroubleshootStep(
                        code="locate_export_sources",
                        description="Identify which buildings inject power during backflow timesteps.",
                        data_needed=["per-building P time series"],
                    ),
                ],
                sentence="Backflow indicates net export to the supply; likely high DER generation during low demand.",
                notes={"backflow_time_fraction": frac},
            )
        )

    # --- Line overload (overcurrent) ---
    if oc:
        frac = _coverage_fraction(oc, T)
        conf = min(1.0, 0.5 + 1.5 * frac + 0.05 * oc_lines)
        evs = [_mk_evidence_from_event(e, "overcurrent") for e in oc]
        findings.append(
            DiagnosisFinding(
                code=DiagnosisCode.LINE_OVERLOAD,
                title="Line overcurrent events detected",
                confidence=conf,
                severity=_severity_from_conf(conf),
                evidence=evs,
                likely_causes=[
                    "cable undersizing vs peak load",
                    "unexpected power routing / topology bottleneck",
                    "export causing reverse loading on feeder",
                ],
                next_steps=[
                    TroubleshootStep(
                        code="rank_lines_by_peak_loading",
                        description="Rank lines by max(I/I_max) and focus on the top offenders.",
                        data_needed=["I time series per line", "line I_max"],
                    ),
                    TroubleshootStep(
                        code="check_topology_bottlenecks",
                        description="Check whether certain branches carry most load due to radial constraints.",
                        data_needed=["network topology", "downstream load allocation"],
                    ),
                ],
                sentence="Some lines exceed their current limits; likely localized congestion or undersized cables.",
                notes={"overcurrent_time_fraction": frac, "lines_with_overcurrent": oc_lines},
            )
        )

    # --- Widespread undervoltage / overvoltage heuristics ---
    # We use counts + time coverage as first-order indicators (v1).
    uv_time = _coverage_fraction(uv, T)
    ov_time = _coverage_fraction(ov, T)

    # “widespread” by node count; if you later store total building node count, this gets better.
    # For now we rely on summary counts only; you can add n_buildings to summary and use fraction.
    # We'll treat >=8 affected nodes as "widespread" if we don't know totals.
    widespread_uv = (uv_nodes >= 8) or (uv_time >= widespread_time_fraction)
    widespread_ov = (ov_nodes >= 8) or (ov_time >= widespread_time_fraction)

    if widespread_uv:
        conf = min(1.0, 0.5 + 3.0 * uv_time + 0.03 * uv_nodes)
        evs = [_mk_evidence_from_event(e, "undervoltage") for e in uv]
        findings.append(
            DiagnosisFinding(
                code=DiagnosisCode.WIDESPREAD_UNDERVOLTAGE,
                title="Widespread undervoltage pattern",
                confidence=conf,
                severity=_severity_from_conf(conf),
                evidence=evs,
                likely_causes=[
                    "upstream voltage setpoint too low / tap position",
                    "insufficient capacity or high R/X drop under peak demand",
                    "voltage regulation absent/ineffective",
                ],
                next_steps=[
                    TroubleshootStep(
                        code="check_transformer_tap_setpoint",
                        description="Check transformer/tap changer setpoints and voltage regulation strategy.",
                        data_needed=["tap changer data", "substation voltage setpoint", "regulator control logic"],
                    ),
                    TroubleshootStep(
                        code="compare_with_peak_demand",
                        description="Check if undervoltage aligns with peak demand periods.",
                        data_needed=["demand profile", "time-of-day labels"],
                    ),
                ],
                sentence="Undervoltage appears widespread/sustained; investigate upstream setpoint and voltage regulation.",
                notes={"uv_nodes": uv_nodes, "uv_time_fraction": uv_time},
            )
        )

    if widespread_ov:
        conf = min(1.0, 0.5 + 3.0 * ov_time + 0.03 * ov_nodes)
        evs = [_mk_evidence_from_event(e, "overvoltage") for e in ov]
        findings.append(
            DiagnosisFinding(
                code=DiagnosisCode.WIDESPREAD_OVERVOLTAGE,
                title="Widespread overvoltage pattern",
                confidence=conf,
                severity=_severity_from_conf(conf),
                evidence=evs,
                likely_causes=[
                    "upstream voltage setpoint too high / tap position",
                    "reactive power control not absorbing VARs",
                    "DER generation raises local voltages",
                ],
                next_steps=[
                    TroubleshootStep(
                        code="check_voltage_setpoint_high",
                        description="Verify substation/tap setpoint is not too high for the feeder conditions.",
                        data_needed=["substation setpoint", "tap settings"],
                    ),
                    TroubleshootStep(
                        code="inspect_DER_voltage_control",
                        description="Check inverter Volt/VAR and Volt/Watt control settings (if applicable).",
                        data_needed=["inverter settings", "DER locations"],
                    ),
                ],
                sentence="Overvoltage appears widespread/sustained; investigate setpoints and DER voltage control.",
                notes={"ov_nodes": ov_nodes, "ov_time_fraction": ov_time},
            )
        )

    # --- Local undervoltage coincident with overcurrent (classic “overload voltage drop”) ---
    # We look for same timesteps between any undervoltage node event and any overcurrent line event.
    if uv and oc:
        best_overlap: Tuple[int, int, List[int]] = (-1, -1, [])
        for uev in uv:
            for cev in oc:
                inter = _intersect_timesteps(uev.timesteps, cev.timesteps)
                if len(inter) > len(best_overlap[2]):
                    best_overlap = (int(uev.where), int(cev.where), inter)



    # --- Local overvoltage with backflow (DER export signature) ---
    if ov and bf:
        bf_ts = _all_timesteps(bf).tolist()
        best_node = (-1, [])
        for oev in ov:
            inter = _intersect_timesteps(oev.timesteps, bf_ts)
            if len(inter) > len(best_node[1]):
                best_node = (int(oev.where), inter)

        if len(best_node[1]) >= coincidence_min_steps:
            node_i, inter = best_node
            frac = (len(inter) / T) if (T and T > 0) else 0.0
            conf = min(1.0, 0.55 + 8.0 * frac)
            findings.append(
                DiagnosisFinding(
                    code=DiagnosisCode.LOCAL_OVERVOLTAGE_WITH_BACKFLOW,
                    title="Overvoltage coincides with reverse power flow",
                    confidence=conf,
                    severity=_severity_from_conf(conf),
                    evidence=[
                        _mk_evidence_from_event(next(e for e in ov if int(e.where) == node_i), "overvoltage"),
                        Evidence(
                            fault_type="backflow",
                            where=-1,
                            timesteps=bf_ts,
                            details={},
                        ),
                        Evidence(
                            fault_type="coincidence",
                            where=-1,
                            timesteps=inter,
                            details={"node": node_i, "overlap_steps": len(inter)},
                        ),
                    ],
                    likely_causes=[
                        "distributed generation export raises local voltage",
                        "inverter voltage control not limiting rise",
                    ],
                    next_steps=[
                        TroubleshootStep(
                            code="identify_exporting_buildings",
                            description="During overlap timesteps, identify top exporting buildings near that node.",
                            data_needed=["per-building export P", "electrical distance / feeder grouping"],
                        ),
                        TroubleshootStep(
                            code="apply_inverter_controls",
                            description="Evaluate Volt/VAR, Volt/Watt, export limiting, or local storage as mitigations.",
                            data_needed=["inverter capability curves", "control settings"],
                        ),
                    ],
                    sentence=f"Overvoltage at node {node_i} overlaps with backflow; likely DER export-driven voltage rise.",
                    notes={"overlap_fraction": frac, "node": node_i},
                )
            )

    # --- Mixed voltage issues: both OV and UV appear (common in grids with DER + peaks) ---
    if ov and uv:
        conf = min(1.0, 0.35 + 2.0 * (ov_time + uv_time))
        findings.append(
            DiagnosisFinding(
                code=DiagnosisCode.MIXED_VOLTAGE_ISSUES,
                title="Both overvoltage and undervoltage occur",
                confidence=conf,
                severity=_severity_from_conf(conf),
                evidence = (
                    [_mk_evidence_from_event(e, "overvoltage") for e in ov[:3]]
                    + [_mk_evidence_from_event(e, "undervoltage") for e in uv[:3]]
                ),
                likely_causes=[
                    "high demand causes voltage drop at peak times",
                    "DER export causes voltage rise at other times",
                    "regulation/control not tuned for both regimes",
                ],
                next_steps=[
                    TroubleshootStep(
                        code="separate_by_time_regime",
                        description="Split analysis into peak-demand windows vs high-DER windows and diagnose each regime separately.",
                        data_needed=["time-of-day labels", "PV/export profiles"],
                    ),
                    TroubleshootStep(
                        code="consider_bidirectional_voltage_regulation",
                        description="Consider coordinated voltage control (tap + inverter VAR + export limits).",
                        data_needed=["regulator models", "inverter control capabilities"],
                    ),
                ],
                sentence="Both OV and UV exist; likely different operating regimes (peak demand vs DER export).",
                notes={"ov_time_fraction": ov_time, "uv_time_fraction": uv_time},
            )
        )

    if not findings:
        findings.append(
            DiagnosisFinding(
                code=DiagnosisCode.UNKNOWN_PATTERN,
                title="No diagnosable pattern from current fault set",
                confidence=0.2,
                severity="info",
                evidence=[],
                likely_causes=["insufficient observables", "thresholds too loose", "faults are rare/noisy"],
                next_steps=[
                    TroubleshootStep(
                        code="add_observables",
                        description="Add more observables (per-node P/Q, voltage at supply, per-line P/Q) to improve diagnosis.",
                        data_needed=["per-node injection", "per-line power flow"],
                    )
                ],
                sentence="Not enough structured evidence to propose a diagnosis confidently.",
            )
        )

    # Report summary: machine-friendly overview
    rep = DiagnosisReport(
        grid_index=grid_index,
        findings=sorted(findings, key=lambda x: (-x.confidence, x.code.value)),
        summary={
            "T": T,
            "fault_types_present": gf.summary.get("fault_types_present", []),
            "num_findings": len(findings),
        },
    )
    return rep


def apply_diagnosis_actions_to_system(
    grid,
    diagnosis_report,
    max_actions: int = 10,
):
    actions = generate_actions_from_diagnosis(grid, diagnosis_report)

    actions = rank_actions(actions)
    actions = actions[:max_actions]

    for action in actions:
        apply_action(grid, action)

    return grid, actions

from dataclasses import dataclass
from typing import Any, Dict, Optional

@dataclass
class Action:
    code: str
    target_type: str          # "grid" | "node" | "line"
    target_id: int            # -1 for grid
    params: Dict[str, Any]    # {"delta": ..., etc}
    score: float = 0.0
    

def generate_actions_from_diagnosis(grid, dr):

    actions = []

    for f in dr.findings:

        if f.code.value == "reverse_power_flow":
            actions += _actions_backflow(grid, f)

        elif f.code.value == "local_overvoltage_with_backflow":
            actions += _actions_ov_backflow(grid, f)

        elif f.code.value == "line_overload":
            actions += _actions_overcurrent(grid, f)

        elif f.code.value == "widespread_undervoltage":
            actions += _actions_widespread_uv(grid, f)

        elif f.code.value == "widespread_overvoltage":
            actions += _actions_widespread_ov(grid, f)

        elif f.code.value == "supply_capacity_exceeded":
            actions += _actions_overflow(grid, f)

    return actions

def _actions_backflow(grid, finding):

    ts = _extract_fault_timesteps(finding)

    exporters = _find_exporting_nodes(grid, ts)

    actions = []

    for node_i, strength in exporters[:5]:
        actions.append(Action(
            code="limit_pv_export",
            target_type="node",
            target_id=node_i,
            params={"factor": 0.7},
            score=strength,
        ))

    return actions

def _actions_ov_backflow(grid, finding):

    ts = _extract_overlap_timesteps(finding)

    exporters = _find_exporting_nodes(grid, ts)

    actions = []

    for node_i, strength in exporters[:5]:
        actions.append(Action(
            code="enable_volt_watt_control",
            target_type="node",
            target_id=node_i,
            params={"slope": 0.05},
            score=strength,
        ))

    return actions

def _actions_overcurrent(grid, finding):

    actions = []

    for ev in finding.evidence:
        if ev.fault_type != "overcurrent":
            continue

        actions.append(Action(
            code="upgrade_line_capacity",
            target_type="line",
            target_id=ev.where,
            params={"factor": 1.3},
            score=1.0,
        ))

    return actions

def _actions_widespread_uv(grid, finding):

    return [
        Action(
            code="increase_slack_voltage",
            target_type="grid",
            target_id=-1,
            params={"delta": 0.02},
            score=finding.confidence,
        )
    ]

def _actions_widespread_ov(grid, finding):

    return [
        Action(
            code="decrease_slack_voltage",
            target_type="grid",
            target_id=-1,
            params={"delta": -0.02},
            score=finding.confidence,
        )
    ]

def _actions_overflow(grid, finding):

    return [
        Action(
            code="increase_supply_capacity",
            target_type="grid",
            target_id=-1,
            params={"factor": 1.2},
            score=finding.confidence,
        )
    ]

def _extract_fault_timesteps(finding):
    ts = []
    for ev in finding.evidence:
        ts.extend(ev.timesteps)
    return sorted(set(ts))


def _extract_overlap_timesteps(finding):
    for ev in finding.evidence:
        if ev.fault_type == "coincidence":
            return ev.timesteps
    return []

def _find_exporting_nodes(grid, timesteps):

    exporters = []

    for i, node in enumerate(grid.nodes):
        P = _to_1d(node.P_inj)
        if P.size == 0:
            continue

        P_fault = P[timesteps]
        export = np.nanmean(np.maximum(P_fault, 0))

        if export > 0:
            exporters.append((i, export))

    exporters.sort(key=lambda x: -x[1])
    return exporters

def rank_actions(actions):
    return sorted(actions, key=lambda a: -a.score)

def apply_action(grid, action):

    if action.code == "increase_slack_voltage":
        _apply_increase_slack_voltage(grid, action)

    elif action.code == "decrease_slack_voltage":
        _apply_decrease_slack_voltage(grid, action)

    elif action.code == "limit_pv_export":
        _apply_limit_pv_export(grid, action)

    elif action.code == "upgrade_line_capacity":
        _apply_upgrade_line_capacity(grid, action)

    elif action.code == "enable_volt_watt_control":
        _apply_volt_watt(grid, action)

    elif action.code == "increase_supply_capacity":
        _apply_increase_supply_capacity(grid, action)
        
def _apply_increase_slack_voltage(grid, action):
    grid.slack_voltage += action.params["delta"]
    
def _apply_decrease_slack_voltage(grid, action):
    grid.slack_voltage -= action.params["delta"]

def _apply_limit_pv_export(grid, action):
    node = grid.nodes[action.target_id]
    node.P_inj = node.P_inj * action.params["factor"]


def _apply_upgrade_line_capacity(grid, action):
    line = grid.lines[action.target_id]
    line.max_current_a *= action.params["factor"]


def _apply_volt_watt(grid, action):
    node = grid.nodes[action.target_id]
    node.volt_watt_enabled = True


def _apply_increase_supply_capacity(grid, action):
    grid.max_supply_apparent_power *= action.params["factor"]  

