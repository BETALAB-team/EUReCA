from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple
import copy
import numpy as np


class FixType(str, Enum):
    SLACK_SET_VOLTAGE = "slack_set_voltage"
    LINE_UPSIZE_AMPACITY = "line_upsize_ampacity"
    LINE_UPSIZE_IMPEDANCE = "line_upsize_impedance"
    CURTAIL_PV = "curtail_pv"
    LIMIT_EXPORT = "limit_export"
    SET_LOAD_PF = "set_load_pf"
    SET_PV_PF = "set_pv_pf"


@dataclass
class FixAction:
    """
    One atomic change to apply to a GridSystem.

    Notes
    -----
    - node actions use node_bus_id (bus_id in Node)
    - line actions use line_id (line_id in Line)
    - value meaning depends on fix_type
    """
    fix_type: FixType
    grid_index: Optional[int] = None
    node_bus_id: Optional[int] = None
    line_id: Optional[int] = None
    value: Any = None
    meta: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TroubleshootPlan:
    """
    A set of FixActions that can be applied to a GridSystem to produce gs_fix.
    """
    title: str
    actions: List[FixAction] = field(default_factory=list)
    tags: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TroubleshootResult:
    """
    Output of running a TroubleshootPlan.

    Fields
    ------
    gs_fix : GridSystem
        Deep-copied and modified system.
    applied : list[FixAction]
        Actions that were successfully applied.
    skipped : list[(FixAction, str)]
        Actions skipped with reason.
    """
    gs_fix: Any
    applied: List[FixAction]
    skipped: List[Tuple[FixAction, str]]


def _iter_grids(gs) -> List[Any]:
    grids = getattr(gs, "grids", None)
    if grids is None:
        grids = gs.split_into_grids()
        gs.grids = grids
    return list(grids)


def apply_troubleshoot_plan(gs, plan: TroubleshootPlan) -> TroubleshootResult:
    """
    Returns a deep-copied GridSystem (gs_fix) with all feasible actions applied.
    """
    gs_fix = copy.deepcopy(gs)
    grids = _iter_grids(gs_fix)

    applied: List[FixAction] = []
    skipped: List[Tuple[FixAction, str]] = []

    for act in plan.actions:
        try:
            gi = act.grid_index
            target_grids = grids if gi is None else [grids[int(gi)]]

            if act.fix_type == FixType.SLACK_SET_VOLTAGE:
                v = float(act.value)
                if not (np.isfinite(v) and v > 0):
                    raise ValueError("value must be a positive finite voltage")
                for g in target_grids:
                    for n in g.nodes:
                        if getattr(n, "node_type", "") == "supply":
                            object.__setattr__(n, "set_voltage", v)
                applied.append(act)
                continue

            if act.fix_type == FixType.LINE_UPSIZE_AMPACITY:
                factor = float(act.value)
                if not (np.isfinite(factor) and factor > 0):
                    raise ValueError("value must be a positive factor")
                for g in target_grids:
                    for l in g.lines:
                        if act.line_id is not None and int(getattr(l, "line_id", 0)) != int(act.line_id):
                            continue
                        mc = getattr(l, "max_current_a", None)
                        if mc is None or not np.isfinite(float(mc)) or float(mc) <= 0:
                            continue
                        object.__setattr__(l, "max_current_a", float(mc) * factor)
                applied.append(act)
                continue

            if act.fix_type == FixType.LINE_UPSIZE_IMPEDANCE:
                factor = float(act.value)
                if not (np.isfinite(factor) and factor > 0):
                    raise ValueError("value must be a positive factor")
                for g in target_grids:
                    for l in g.lines:
                        if act.line_id is not None and int(getattr(l, "line_id", 0)) != int(act.line_id):
                            continue
                        r = getattr(l, "r_ohm_per_km", None)
                        x = getattr(l, "x_ohm_per_km", None)
                        if r is not None and np.isfinite(float(r)) and float(r) > 0:
                            object.__setattr__(l, "r_ohm_per_km", float(r) / factor)
                        if x is not None and np.isfinite(float(x)) and float(x) > 0:
                            object.__setattr__(l, "x_ohm_per_km", float(x) / factor)
                applied.append(act)
                continue

            if act.fix_type == FixType.CURTAIL_PV:
                alpha = float(act.value)
                if not (np.isfinite(alpha) and 0 <= alpha <= 1):
                    raise ValueError("value must be in [0, 1]")
                for g in target_grids:
                    for n in g.nodes:
                        if getattr(n, "node_type", "") != "building":
                            continue
                        if act.node_bus_id is not None and int(getattr(n, "bus_id", 0)) != int(act.node_bus_id):
                            continue
                        bp = getattr(n, "building_production", None)
                        if bp is None:
                            continue
                        bp2 = np.asarray(bp, dtype=float) * alpha
                        object.__setattr__(n, "building_production", bp2)
                applied.append(act)
                continue

            if act.fix_type == FixType.LIMIT_EXPORT:
                pmax_w = float(act.value)
                if not (np.isfinite(pmax_w) and pmax_w >= 0):
                    raise ValueError("value must be >= 0 (W)")
                for g in target_grids:
                    for n in g.nodes:
                        if getattr(n, "node_type", "") != "building":
                            continue
                        if act.node_bus_id is not None and int(getattr(n, "bus_id", 0)) != int(act.node_bus_id):
                            continue
                        bp = getattr(n, "building_production", None)
                        if bp is None:
                            continue
                        bc = getattr(n, "building_consumption", None)
                        bc = np.zeros_like(bp, dtype=float) if bc is None else np.asarray(bc, dtype=float).reshape(-1)
                        bp = np.asarray(bp, dtype=float).reshape(-1)
                        T = min(bp.size, bc.size)
                        bp = bp[:T]
                        bc = bc[:T]
                        inj = np.maximum(bp - bc, 0.0)
                        inj_limited = np.minimum(inj, pmax_w)
                        bp2 = inj_limited + bc
                        object.__setattr__(n, "building_production", bp2)
                applied.append(act)
                continue

            if act.fix_type == FixType.SET_LOAD_PF:
                pf = float(act.value)
                if not (np.isfinite(pf) and 0 < pf <= 1):
                    raise ValueError("pf must be in (0, 1]")
                for g in target_grids:
                    for n in g.nodes:
                        if getattr(n, "node_type", "") != "building":
                            continue
                        if act.node_bus_id is not None and int(getattr(n, "bus_id", 0)) != int(act.node_bus_id):
                            continue
                        object.__setattr__(n, "pf_load", pf)
                applied.append(act)
                continue

            if act.fix_type == FixType.SET_PV_PF:
                pf = float(act.value)
                if not (np.isfinite(pf) and 0 < pf <= 1):
                    raise ValueError("pf must be in (0, 1]")
                for g in target_grids:
                    for n in g.nodes:
                        if getattr(n, "node_type", "") != "building":
                            continue
                        if act.node_bus_id is not None and int(getattr(n, "bus_id", 0)) != int(act.node_bus_id):
                            continue
                        object.__setattr__(n, "pf_pv", pf)
                applied.append(act)
                continue

            skipped.append((act, f"unknown fix_type={act.fix_type!s}"))
        except Exception as e:
            skipped.append((act, repr(e)))

    return TroubleshootResult(gs_fix=gs_fix, applied=applied, skipped=skipped)


def _get(obj, *names, default=None):
    for n in names:
        if hasattr(obj, n):
            return getattr(obj, n)
    if isinstance(obj, dict):
        for n in names:
            if n in obj:
                return obj[n]
    return default


def _severity_weight(sev: str) -> float:
    s = (sev or "").lower()
    if s == "critical":
        return 3.0
    if s == "warning":
        return 2.0
    return 1.0


def _extract_events_from_diagnosis(diagnosis) -> List[Dict[str, Any]]:
    """
    Returns list of normalized event dicts:
      {
        "code": str,
        "severity": str,
        "nodes": [(bus_id, timesteps_count, timesteps_list)],
        "lines": [(line_id, timesteps_count, timesteps_list)],
        "supply": [(key, timesteps_count, timesteps_list)]
      }
    """

    findings = _get(diagnosis, "findings", default=[])
    out: List[Dict[str, Any]] = []

    for f in findings:
        code = str(_get(f, "code", "name", default="UNKNOWN"))
        sev = str(_get(f, "severity", default="info"))
        ev_list = list(_get(f, "evidence", default=[]))

        nodes: List[Tuple[int, int, List[int]]] = []
        lines: List[Tuple[int, int, List[int]]] = []
        supply: List[Tuple[str, int, List[int]]] = []

        for e in ev_list:
            
            kind = str(_get(e, "fault_type", "type", default="")).lower()
            
            ts = _get(e, "timesteps", "steps", default=[])

            ts = [] if ts is None else list(ts)
            tcnt = len(ts)
            if "current" in kind:
                bus_id = None
                line_id = _get(e, "where", default=None)
                
            if "voltage" in kind:
                line_id = None
                bus_id = _get(e, "where", default=None)
            
            if "flow" in kind:
                line_id = None
                bus_id = _get(e, "where", default=None)
                
            if bus_id is not None or "node" in kind:
                if bus_id is None:
                    bus_id = _get(e, "id", default=None)
                if bus_id is not None:
                    nodes.append((int(bus_id), tcnt, ts))
                    continue

            if line_id is not None or "line" in kind:
                if line_id is None:
                    line_id = _get(e, "id", default=None)
                if line_id is not None:
                    lines.append((int(line_id), tcnt, ts))
                    continue

            if "supply" in kind or code.lower().startswith("supply"):
                supply.append(("supply", tcnt, ts))

        out.append({"code": code, "severity": sev, "nodes": nodes, "lines": lines, "supply": supply})
    return out


def _rank_by_timesteps(items: List[Tuple[int, int, List[int]]], top_k: int = 3) -> List[Tuple[int, int, List[int]]]:
    return sorted(items, key=lambda x: x[1], reverse=True)[:top_k]


def _plan_from_events(
    grid,
    grid_index: int,
    events: List[Dict[str, Any]],
    *,
    min_voltage: float,
    max_voltage: float,
    max_supply_S: float,
    tap_step_v: float = 2.0,
    tap_max_steps: int = 4,
    pv_curtail_alpha: float = 0.85,
    export_limit_w: float = 2000.0,
    line_ampacity_factor: float = 1.25,
    line_impedance_factor: float = 1.25,
) -> TroubleshootPlan:
    """
    Build a plan for ONE grid based on diagnosis events.
    Always returns a TroubleshootPlan (possibly with zero actions).
    """
    actions: List[FixAction] = []
    meta: Dict[str, Any] = {"grid_index": grid_index, "rules": {}, "events": [e.get("code") for e in events]}

    slack_v = None
    for n in getattr(grid, "nodes", []):
        if getattr(n, "node_type", "") == "supply":
            sv = getattr(n, "set_voltage", None)
            if sv is not None and np.isfinite(float(sv)) and float(sv) > 0:
                slack_v = float(sv)
                break
    if slack_v is None:
        slack_v = float(max_voltage)

    need_tap_up = 0.0
    need_tap_down = 0.0
    overloaded_lines: List[Tuple[int, int, List[int], float]] = []
    ov_nodes: List[Tuple[int, int, List[int], float]] = []
    uv_nodes: List[Tuple[int, int, List[int], float]] = []
    backflow_nodes: List[Tuple[int, int, List[int], float]] = []

    for ev in events:
        code = str(ev.get("code", "")).lower()
        w = _severity_weight(str(ev.get("severity", "info")))

        top_lines = _rank_by_timesteps(ev.get("lines", []), top_k=5)
        top_nodes = _rank_by_timesteps(ev.get("nodes", []), top_k=5)

        if "undervoltage" in code or "low_voltage" in code:
            need_tap_up += 1.0 * w
            for bus_id, cnt, ts in top_nodes:
                uv_nodes.append((bus_id, cnt, ts, w))

        if "overvoltage" in code or "high_voltage" in code:
            need_tap_down += 1.0 * w
            for bus_id, cnt, ts in top_nodes:
                ov_nodes.append((bus_id, cnt, ts, w))

        if "overcurrent" in code or "congestion" in code or "line_overload" in code:
            for lid, cnt, ts in top_lines:
                overloaded_lines.append((lid, cnt, ts, w))

        if "backflow" in code or "reverse_power" in code or "export" in code:
            for bus_id, cnt, ts in top_nodes:
                backflow_nodes.append((bus_id, cnt, ts, w))

        if "overflow" in code or "supply_overflow" in code or "supply_limit" in code:
            meta["rules"]["supply_overflow_seen"] = True

    tap_score = need_tap_up - need_tap_down
    if abs(tap_score) > 0.5:
        steps = int(min(tap_max_steps, max(1, round(abs(tap_score)))))
        dv = tap_step_v * steps * (1.0 if tap_score > 0 else -1.0)
        new_slack_v = float(np.clip(slack_v + dv, min_voltage, max_voltage))
        actions.append(
            FixAction(
                fix_type=FixType.SLACK_SET_VOLTAGE,
                grid_index=grid_index,
                value=new_slack_v,
                meta={"from": slack_v, "dv": dv, "score": tap_score},
            )
        )

    if overloaded_lines:
        ranked = sorted(overloaded_lines, key=lambda x: x[1] * x[3], reverse=True)[:3]
        for lid, cnt, ts, w in ranked:
            actions.append(
                FixAction(
                    fix_type=FixType.LINE_UPSIZE_AMPACITY,
                    grid_index=grid_index,
                    line_id=lid,
                    value=line_ampacity_factor,
                    meta={"hits": cnt, "severity_weight": w},
                )
            )
            actions.append(
                FixAction(
                    fix_type=FixType.LINE_UPSIZE_IMPEDANCE,
                    grid_index=grid_index,
                    line_id=lid,
                    value=line_impedance_factor,
                    meta={"hits": cnt, "severity_weight": w},
                )
            )

    if ov_nodes or backflow_nodes:
        merged: Dict[int, Dict[str, Any]] = {}
        for bus_id, cnt, ts, w in (ov_nodes + backflow_nodes):
            s = merged.get(int(bus_id), {"cnt": 0, "w": 0.0, "ts": []})
            s["cnt"] += int(cnt)
            s["w"] += float(w)
            s["ts"] = (s["ts"] + list(ts))[:50]
            merged[int(bus_id)] = s

        ranked_nodes = sorted(merged.items(), key=lambda kv: kv[1]["cnt"] * kv[1]["w"], reverse=True)[:5]
        for bus_id, s in ranked_nodes:
            actions.append(
                FixAction(
                    fix_type=FixType.LIMIT_EXPORT,
                    grid_index=grid_index,
                    node_bus_id=int(bus_id),
                    value=float(export_limit_w),
                    meta={"hits": s["cnt"], "severity_weight": s["w"]},
                )
            )
            actions.append(
                FixAction(
                    fix_type=FixType.CURTAIL_PV,
                    grid_index=grid_index,
                    node_bus_id=int(bus_id),
                    value=float(pv_curtail_alpha),
                    meta={"hits": s["cnt"], "severity_weight": s["w"]},
                )
            )

    meta["rules"]["max_supply_S_seen"] = float(max_supply_S)

    return TroubleshootPlan(
        title=f"Auto plan for Grid {grid_index}",
        actions=actions,
        tags=meta,
    )


def build_troubleshoot_plan_auto(
    gs,
    diagnoses_by_grid: Dict[int, Any],
    *,
    min_voltage: float,
    max_voltage: float,
    max_supply_apparent_power: float,
    tap_step_v: float = 2.0,
    tap_max_steps: int = 4,
    pv_curtail_alpha: float = 0.85,
    export_limit_w: float = 2000.0,
    line_ampacity_factor: float = 1.25,
    line_impedance_factor: float = 1.25,
) -> TroubleshootPlan:
    grids = _iter_grids(gs)

    all_actions: List[FixAction] = []
    tags: Dict[str, Any] = {"auto": True, "per_grid": {}}

    for gi, grid in enumerate(grids):
        diag_obj = diagnoses_by_grid.get(gi, None)
        if diag_obj is None:
            tags["per_grid"][gi] = {"n_actions": 0, "events": [], "note": "no diagnosis"}
            continue

        events = _extract_events_from_diagnosis(diag_obj)
        plan_g = _plan_from_events(
            grid,
            gi,
            events,
            min_voltage=min_voltage,
            max_voltage=max_voltage,
            max_supply_S=max_supply_apparent_power,
            tap_step_v=tap_step_v,
            tap_max_steps=tap_max_steps,
            pv_curtail_alpha=pv_curtail_alpha,
            export_limit_w=export_limit_w,
            line_ampacity_factor=line_ampacity_factor,
            line_impedance_factor=line_impedance_factor,
        )

        all_actions.extend(plan_g.actions)
        tags["per_grid"][gi] = {"n_actions": len(plan_g.actions), "events": [e.get("code") for e in events]}

    return TroubleshootPlan(
        title="Auto troubleshoot plan (all grids)",
        actions=all_actions,
        tags=tags,
    )


def run_troubleshoot(
    gs,
    plan: TroubleshootPlan,
    *,
    analyze_topology_fn,
    apply_results_to_grid_fn,
    fault_report_fn,
    diagnose_fn,
    solve_kwargs: Optional[Dict[str, Any]] = None,
):
    solve_kwargs = {} if solve_kwargs is None else dict(solve_kwargs)

    tr = apply_troubleshoot_plan(gs, plan)
    gs_fix = tr.gs_fix

    grids = _iter_grids(gs_fix)
    for g in grids:
        g.kind = analyze_topology_fn(g)
        g.solve(**solve_kwargs)
        apply_results_to_grid_fn(g, g.results)

    gf_fix = fault_report_fn(gs_fix)
    diagnoses_fix = {i: diagnose_fn(fault, i) for i, fault in gf_fix.grids.items()}

    return tr, gf_fix, diagnoses_fix
