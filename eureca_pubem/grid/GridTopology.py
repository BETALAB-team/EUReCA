from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import numpy as np
import warnings


@dataclass
class TopologyInfo:
    is_connected: bool
    N: int
    E: int
    loops: int
    kind: str
    slack_bus: int
    drop_line_ids: List[int]
    tree_line_ids: List[int]


def _pick_slack_bus(grid, slack_bus_id):
    supply = [int(n.bus_id) for n in grid.nodes if n.node_type == "supply"]
    if slack_bus_id is not None and int(slack_bus_id) in supply:
        return int(slack_bus_id)
    if supply:
        return int(supply[0])
    return int(grid.nodes[0].bus_id)


def _build_adjacency(lines):
    adj = {}
    for l in lines:
        u = int(l.from_bus)
        v = int(l.to_bus)
        lid = int(l.line_id)
        adj.setdefault(u, []).append((v, lid))
        adj.setdefault(v, []).append((u, lid))
    return adj


def _spanning_tree_lines(grid, slack_bus):
    buses = [int(n.bus_id) for n in grid.nodes]
    bus_set = set(buses)
    adj = _build_adjacency(grid.lines)

    parent = {slack_bus: None}
    parent_line = {slack_bus: None}
    q = [slack_bus]
    seen = {slack_bus}

    while q:
        u = q.pop(0)
        for v, lid in adj.get(u, []):
            if v not in bus_set:
                continue
            if v in seen:
                continue
            seen.add(v)
            parent[v] = u
            parent_line[v] = lid
            q.append(v)

    is_connected = (len(seen) == len(bus_set))
    tree_line_ids = [int(parent_line[b]) for b in seen if b != slack_bus and parent_line.get(b) is not None]

    return is_connected, tree_line_ids


def analyze_topology(
    grid,
    *,
    slack_bus_id: Optional[int] = None,
    weakly_meshed_max_loops: int = 2,
    weakly_meshed_max_ratio: float = 0.02,
) -> TopologyInfo:
    slack_bus = _pick_slack_bus(grid, slack_bus_id)
    N = len(grid.nodes)
    E = len(grid.lines)
    loops = int(E - (N - 1))

    is_connected, tree_line_ids = _spanning_tree_lines(grid, slack_bus)
    if not is_connected:
        warnings.warn("Grid not fully connected from slack, classification may be unreliable")

    tree_line_set = set(tree_line_ids)
    all_line_ids = [int(l.line_id) for l in grid.lines]
    drop_line_ids = [lid for lid in all_line_ids if lid not in tree_line_set]

    ratio = (loops / max(N, 1)) if loops > 0 else 0.0
    if (not is_connected) and N > 0:
        kind = "disconnected"
    elif loops <= 0:
        kind = "radial"
        drop_line_ids = []
    elif loops <= weakly_meshed_max_loops or ratio <= weakly_meshed_max_ratio:
        kind = "weakly_meshed"
    else:
        kind = "meshed"

    return TopologyInfo(
        is_connected=is_connected,
        N=N,
        E=E,
        loops=max(loops, 0),
        kind=kind,
        slack_bus=slack_bus,
        drop_line_ids=drop_line_ids,
        tree_line_ids=tree_line_ids,
    )
