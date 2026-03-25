from dataclasses import dataclass, field
from typing import Dict, List, Any, Union


@dataclass
class Building:
    building_id: int

    gas_consumed: float = 0.0
    energy_bought_from_dhn: float = 0.0

    electricity_bought: float = 0.0
    electricity_sold: float = 0.0

    voltage: Union[float, None] = None
    dhn_supply_pressure: Union[float, None] = None
    dhn_supply_temperature: Union[float, None] = None



@dataclass
class ElectricalGrid:
    grid_id: int
    nodes: Dict[int, Any] = field(default_factory=dict)
    lines: Dict[int, Any] = field(default_factory=dict)
    results: Union[float, None] = None


@dataclass
class ElectricalDiagnosis:
    grid_id: int
    faults: List[Any] = field(default_factory=list)



@dataclass
class DHNSubsystem:
    subsystem_id: int
    nodes: Dict[int, Any] = field(default_factory=dict)
    pipes: Dict[int, Any] = field(default_factory=dict)
    results: Union[float, None] = None


@dataclass
class DHNDiagnosis:
    subsystem_id: int
    faults: List[Any] = field(default_factory=list)



@dataclass
class Scenario:
    buildings: Dict[int, Building] 
    Electrical_Network: dict
    District_Heating_Systems: list
    Electrical_Diagnoses: list
    District_Heating_Diagnoses: list


