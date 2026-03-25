from __future__ import annotations

from dataclasses import dataclass 
from typing import  List
import numpy as np 


@dataclass 
class Node: 
    node_id: int                    # id to identify the nodes later.
    node_type: str                  # type of node, it can be consumer, supplier, or connection
    demand: np.ndarray   =None           # hourly demand of the node if there is any 
    

@dataclass
class Line:
    line_id: int                          # id to identify the line
    start_node: int                       # id of the node of the start
    end_node: int                         # id of the node of the end
    length: float                         # length of the pipe [m]
    inner_diameter: float     = 1           # inner diameter of the pipe [mm]
    outer_diameter: float   = 1            # outer diameter of the pipe [mm]
    conductivity: float        = 1           # conductivity of the pipe material [W/mK]
    insulation_conductivity: float   = 0.03     # conductivity of the insulation material [W/mK]
    insulation_thickness: float = 5          # thickness of insulation [mm] 
    roughness: float= 0.00001
    pipe_specific_heat_capacity: float= 0.5
    m: float= 0.1
    dp: float= 10
    knees: int = 0
    max_pressure: float | None = None
    
class DistrictHeating:
    def __init__(self, name: str, nodes: List[Node], lines: List[Line]):
        self.name = name
        self.nodes = list(nodes)
        self.lines = list(lines)
        self.active = any(n.node_type == "supply" for n in self.nodes)

