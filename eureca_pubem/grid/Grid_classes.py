from __future__ import annotations

from eureca_pubem.grid.powerflow_timeseries import solve_grid_timeseries_auto

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np

@dataclass
class Inverter:
    """
    Represents a PV inverter attached to a node.

    Parameters
    ----------
    capacity_w : float
        Inverter nominal capacity in W.
    efficiency : float
        Inverter conversion efficiency in [0, 1].
    """
    capacity_w: float = 10_000.0
    efficiency: float = 0.98

    def __post_init__(self):
        if self.capacity_w <= 0:
            raise ValueError("Inverter capacity_w must be > 0")
        if not (0 < self.efficiency <= 1):
            raise ValueError("Inverter efficiency must be in (0, 1]")

@dataclass
class Node:
    """
    Represents a node in the grid system.

    Parameters
    ----------
    bus_id : int
        Unique identifier of the node.
    node_type : str
        Type of node. Allowed values: "supply", "building", "connection".
    x : float
        X coordinate in meters.
    y : float
        Y coordinate in meters.
    building_id : int | None
        Associated building ID if node_type is "building".
    max_power : float | None
        Substation maximum power for supply nodes.
    set_voltage : float | None
        Substation set voltage for supply nodes.
    set_voltage_angle : float | None
        Substation voltage angle for supply nodes.
    building_consumption : np.ndarray | None
        Time series consumption for building nodes.
    building_production : np.ndarray | None
        Time series production for building nodes.
    """
    bus_id: int
    node_type: str
    x: float
    y: float
    building_id: Optional[int] = None
    max_power: Optional[float] = None
    set_voltage: Optional[float] = None
    set_voltage_angle: Optional[float] = None
    building_consumption: Optional[np.ndarray] = None
    building_production: Optional[np.ndarray] = None
    inverter: Optional[Inverter] = None
    pf_load: float = 0.95
    load_lagging: bool = True
    pf_pv: float = 1.0
    pv_lagging: bool = False

    def __post_init__(self):
        t = self.node_type.lower()
        object.__setattr__(self, "node_type", t)

        if t == "building" and self.building_id is None:
            raise ValueError(f"Building node {self.bus_id} requires a building_id")

        if t != "supply" and any(v is not None for v in (self.max_power, self.set_voltage, self.set_voltage_angle)):
            raise ValueError(f"Non-supply node {self.bus_id} cannot have substation parameters")

        if t != "building" and (self.building_consumption is not None or self.building_production is not None):
            raise ValueError(f"Non-building node {self.bus_id} cannot have building consumption/production arrays")
        if self.node_type != "building" and self.inverter is not None:
            raise ValueError(f"Non-building node {self.bus_id} cannot have an inverter")
@dataclass
class Line:
    """
    Represents a line connecting two nodes.

    Parameters
    ----------
    line_id : int
        Unique identifier of the line.
    from_bus : int
        Origin node bus ID.
    to_bus : int
        Destination node bus ID.
    cable_type : int | None
        Cable type identifier (Cable ID from linetypes sheet).
    length_m : float
        Length of the line in meters.
    r_ohm_per_km : float | None
        Resistance in ohm/km for the cable type.
    x_ohm_per_km : float | None
        Reactance in ohm/km for the cable type.
    max_current_a : float | None
        Maximum current in A for the cable type.
    """
    line_id: int
    from_bus: int
    to_bus: int
    length_m: float
    cable_type: Optional[int] = None
    r_ohm_per_km: Optional[float] = None
    x_ohm_per_km: Optional[float] = None
    max_current_a: Optional[float] = None

    origin_line_id: Optional[int] = None
    required_current: Optional[float] = None
    parallel_count: int = 1


@dataclass
class Grid:
    """
    Represents an isolated grid subsystem containing at least one supply node.

    Parameters
    ----------
    nodes : list[Node]
        Nodes belonging to the grid.
    lines : list[Line]
        Lines belonging to the grid.
    """
    nodes: List[Node]
    lines: List[Line]
    def solve(self, timerange=None, **kwargs):
        self.results = solve_grid_timeseries_auto(self, **kwargs)
        
        return self.results
    
class GridSystem:
    def __init__(self, nodes, lines):
        self.nodes = nodes
        self.lines = lines

    def split_into_grids(self):
        adjacency = {bid: set() for bid in self.nodes}
        for l in self.lines:
            adjacency[l.from_bus].add(l.to_bus)
            adjacency[l.to_bus].add(l.from_bus)

        visited = set()
        grids = []

        for start in self.nodes:
            if start in visited:
                continue

            stack = [start]
            component = set()

            while stack:
                n = stack.pop()
                if n in visited:
                    continue
                visited.add(n)
                component.add(n)
                stack.extend(adjacency[n] - visited)

            nodes_comp = [self.nodes[i] for i in component]
            lines_comp = [l for l in self.lines if l.from_bus in component and l.to_bus in component]

            if any(n.node_type == "supply" for n in nodes_comp):
                grids.append(Grid(nodes=nodes_comp, lines=lines_comp))

        return grids