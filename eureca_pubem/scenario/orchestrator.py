from eureca_pubem.dhn.DH_Classes import District_Heating_System
from eureca_pubem.dhn.DH_tools import size_dhn_with_baseline
from eureca_pubem.dhn.DH_diagnosis import diagnose_dhn_faults
from eureca_pubem.dhn.DH_design import attach_missing_dhn_consumers, read_peak_demands, apply_dhn_activation_from_demands
import eureca_pubem.dhn.DH_tools as toolbox
import eureca_pubem.dhn.DH_tools_fast as fast_toolbox

from eureca_pubem.grid.GridGenerator import GridGenerator
from eureca_pubem.grid.GridTopology import analyze_topology
from eureca_pubem.grid.powerflow_timeseries import  apply_results_to_grid, plot_interactive

from time import time
from eureca_pubem.grid import Diagnosis as diag
from eureca_pubem.grid import Troubleshoot as ts
from eureca_pubem.scenario.scenario_objects import Building, Scenario

import numpy as np
import copy
import pandas as pd

def interpolate_1d(arr):
    n = len(arr)
    x = np.arange(n)
    good = ~np.isnan(arr)

    if good.sum() > 1:
        arr = np.interp(x, x[good], arr[good])

    return arr


def interpolate_2d(arr):
    n, m = arr.shape
    x = np.arange(n)

    for col in range(m):
        col_data = arr[:, col]
        good = ~np.isnan(col_data)

        if good.sum() > 1:
            arr[:, col] = np.interp(x, x[good], col_data[good])

    return arr
def compute_bad_timesteps(arr, z_thresh=4.0):
    arr = arr.astype(float)

    median = np.nanmedian(arr)
    mad = np.nanmedian(np.abs(arr - median))

    if mad == 0:
        return np.zeros_like(arr, dtype=bool)

    robust_z = 0.6745 * (arr - median) / mad
    mask = np.abs(robust_z) > z_thresh

    return mask
def clean_with_mask(arr, mask):
    arr = arr.astype(float).copy()

    # if mask is 2D → reduce to 1D timestep mask
    if mask.ndim == 2:
        mask = np.any(mask, axis=1)

    # ensure mask is 1D
    if mask.ndim != 1:
        raise ValueError("Mask must be 1D or 2D")

    # if arr is 1D
    if arr.ndim == 1:
        arr[mask] = np.nan
        return interpolate_1d(arr)

    # if arr is 2D
    elif arr.ndim == 2:
        arr[mask, :] = np.nan
        return interpolate_2d(arr)

    return arr

def compute_bad_timesteps(arr, z_thresh=4.0):
    arr = arr.astype(float)

    median = np.nanmedian(arr)
    mad = np.nanmedian(np.abs(arr - median))

    if mad == 0:
        return np.zeros_like(arr, dtype=bool)

    robust_z = 0.6745 * (arr - median) / mad
    mask = np.abs(robust_z) > z_thresh

    return mask



def clean_object_timeseries(obj, mask):
    """
    Cleans all numpy arrays inside an object that
    have first dimension equal to len(mask).
    """

    T = len(mask)

    for attr, value in obj.__dict__.items():

        if isinstance(value, np.ndarray):

            if value.ndim == 1 and value.shape[0] == T:
                obj.__dict__[attr] = clean_with_mask(value, mask)

            elif value.ndim == 2 and value.shape[0] == T:
                obj.__dict__[attr] = clean_with_mask(value, mask)

        elif isinstance(value, list):
            for v in value:
                if hasattr(v, "__dict__"):
                    clean_object_timeseries(v, mask)

        elif hasattr(value, "__dict__"):
            clean_object_timeseries(value, mask)

def interpolate_inplace(arr):
    n = arr.shape[0]
    if n == 0:
        return

    mask = np.isnan(arr)
    if not mask.any():
        return

    good = ~mask
    if good.sum() < 2:
        return

    x = np.arange(n)
    arr[:] = np.interp(x, x[good], arr[good])
    
    
def clean_scenario_fast(scenario):

    # ---------------- DHN ----------------
    for dhn in scenario.District_Heating_Systems.values():

        ref = getattr(dhn, "hourly_heat_loss", None)
        if ref is None:
            continue

        mask = compute_bad_timesteps(ref)
        T = len(mask)

        # ---- DHN object level ----
        for attr in [
            "hourly_heat_loss",
            "hourly_heat_demand",
        ]:
            arr = getattr(dhn, attr, None)
            if isinstance(arr, np.ndarray) and arr.shape[0] == T:
                arr[mask] = np.nan
                interpolate_inplace(arr)

        # ---- Lines ----
        for line in dhn.lines:
            for attr in [
                "m_init",
                "return_end_temperature",
                "return_loss",
                "return_pressure_drop",
                "return_start_temperature",
                "supply_end_temperature",
                "supply_loss",
                "supply_pressure_drop",
                "supply_start_temperature",
            ]:
                arr = getattr(line, attr, None)
                if isinstance(arr, np.ndarray) and arr.shape[0] == T:
                    arr[mask] = np.nan
                    interpolate_inplace(arr)

        # ---- Nodes ----
        for node in dhn.nodes:
            for attr in [
                "demand",
                "heat_injected",
                "m_init",
                "pump_power",
                "return_pressure",
                "return_temperature",
                "supply_pressure",
                "supply_temperature",
            ]:
                arr = getattr(node, attr, None)
                if isinstance(arr, np.ndarray) and arr.shape[0] == T:
                    arr[mask] = np.nan
                    interpolate_inplace(arr)

    # ---------------- GRID ----------------
    for grid in scenario.Electrical_Network:

        ref = getattr(grid.results, "P_supply", None)
        if ref is None:
            continue

        mask = compute_bad_timesteps(ref)
        mask = np.asarray(mask).reshape(-1)
        T = len(mask)

        # ---- Grid results ----
        for attr in [
            "E_loss_line_kWh",
            "I_line",
            "P_inj",
            "P_loss_line",
            "P_supply",
            "Q_inj",
            "Q_supply",
            "S_supply",
            "V",
            "Va_rad",
            "Vm_pu",
            "converged",
        ]:
            arr = getattr(grid.results, attr, None)
            if isinstance(arr, np.ndarray) and arr.shape[0] == T:
            
                if arr.ndim == 1:
                    arr[mask] = np.nan
                    interpolate_inplace(arr)
            
                elif arr.ndim == 2:
                    arr[mask, :] = np.nan
                    for col in range(arr.shape[1]):
                        interpolate_inplace(arr[:, col])

        # ---- Grid lines ----
        for line in grid.lines:
            for attr in [
                "E_loss_kWh",
                "I",
                "P_loss",
                "current_ts",
            ]:
                arr = getattr(line, attr, None)
                if isinstance(arr, np.ndarray) and arr.shape[0] == T:
                    arr[mask] = np.nan
                    interpolate_inplace(arr)

        # ---- Grid nodes ----
        for node in grid.nodes:
            for attr in [
                "P_inj",
                "Q_inj",
                "V",
                "building_consumption",
                "building_production",
                "consumption",
                "downstream_power",
                "net_power",
                "production",
            ]:
                arr = getattr(node, attr, None)
                if isinstance(arr, np.ndarray) and arr.shape[0] == T:
                    arr[mask] = np.nan
                    interpolate_inplace(arr)

def clean_array(arr, z_thresh=4.0):

    arr = arr.astype(float)

    # ---- 1D ----
    if arr.ndim == 1:
        return _clean_1d(arr, z_thresh)

    # ---- 2D ----
    elif arr.ndim == 2:
        cleaned = np.zeros_like(arr)
        for col in range(arr.shape[1]):
            cleaned[:, col] = _clean_1d(arr[:, col], z_thresh)
        return cleaned

    else:
        return arr


def _clean_1d(arr, z_thresh):

    s = pd.Series(arr)

    # interpolate NaN
    s = s.interpolate(limit_direction="both")

    median = s.median()
    mad = np.median(np.abs(s - median))

    if mad != 0:
        robust_z = 0.6745 * (s - median) / mad
        s[np.abs(robust_z) > z_thresh] = np.nan
        s = s.interpolate(limit_direction="both")

    return s.to_numpy()


def clean_all_arrays(obj):
    if isinstance(obj, np.ndarray):
        return clean_array(obj)

    elif isinstance(obj, dict):
        for k, v in obj.items():
            obj[k] = clean_all_arrays(v)
        return obj

    elif isinstance(obj, (list, tuple)):
        return type(obj)(clean_all_arrays(v) for v in obj)

    elif hasattr(obj, "__dict__"):
        for attr, value in obj.__dict__.items():
            setattr(obj, attr, clean_all_arrays(value))
        return obj

    else:
        return obj
 



    
    
    
def apply_grid_changes(
    scenario,
    consumption_csv_path: str,
    production_csv_path: str,
    baseline_obj,
    *,
    slack_voltage: float = 400.0,
    cable_confidence_factor: float = 1.0,
    timestep_hours: float = 1.0,
):
    """
    Electrical equivalent of apply_dhn_changes.

    Workflow:
        1. Deep-copy scenario
        2. Update electrical profiles
        3. Recompute currents (pre-resize)
        4. Resize against baseline
        5. Solve power flow
        6. Store line changes
    """

    s_2 = scenario

    update_electrical_profiles(
        s_2,
        consumption_csv_path,
        production_csv_path,
    )

    recompute_grid_currents(
        s_2,
        slack_voltage=slack_voltage,
        timestep_hours=timestep_hours,
    )

    resize_grid_against_baseline(
        baseline_obj=baseline_obj,
        scenario=s_2,
        cable_confidence_factor=cable_confidence_factor,
    )

    solve_electrical_network(
        s_2,
        slack_voltage=slack_voltage,
        timestep_hours=timestep_hours,
    )

    return s_2
     
def apply_dhn_changes(scenario, csv_path, baseline_obj, climate_path):
    peak = read_peak_demands(csv_path) 
    s_2 = scenario

    attach_missing_dhn_consumers (s_2, peak, csv_path)

    compare_dhn(baseline_obj, s_2)

    # solve_dhn(s_2.District_Heating_Systems)
    dhs = District_Heating_System.from_subsystems(
    subsystems=s_2.District_Heating_Systems,
    random_seed=17,
    climate_file_path=climate_path
    )
    dhs.solve(time_frame = range(1,8760), #range in hour
                return_pressure = 150_000.0, #Pa
                min_pressure_difference = 100_000, #Pa
                supply_temperature = 20, #C
                required_temperature = 60, #C
                ground_temperature_profile = 7.2, #C
                solve_no_loss = False,
                use_numba = True)
    s_2.District_Heating_Systems = dhs.subsystems
    # s_2.Distsolve(time_frame = range(1,8760), #range in hour
    #             return_pressure = 150_000.0, #Pa
    #             min_pressure_difference = 100_000, #Pa
    #             supply_temperature = 20, #C
    #             required_temperature = 60, #C
    #             ground_temperature_profile = 7.2, #C
    #             solve_no_loss = False,
    #             use_numba = True)
    return s_2

def ensure_origin_line_id(system):
    """
    Ensures every line has origin_line_id.
    Baseline lines get origin_line_id = line_id.
    """
    for l in system.lines:
        if not hasattr(l, "origin_line_id"):
            l.origin_line_id = l.line_id

    
def compare_dhn(baseline_obj, scenario):
    new_pipes = []
    pipes_changed = []
    baseline = baseline_obj
    
    for dhn in baseline.District_Heating_Systems.values():
        ensure_origin_line_id(dhn)
    
    for dhn in scenario.District_Heating_Systems.values():
        ensure_origin_line_id(dhn)
    baseline_lookup = {
        l.origin_line_id: l
        for dhn in baseline.District_Heating_Systems.values()
        for l in dhn.lines
        if l.origin_line_id is not None
    }

    for dhn in scenario.District_Heating_Systems.values():
        for l in dhn.lines:

            if l.origin_line_id is None:
                new_pipes.append(l)
                continue

            base = baseline_lookup.get(l.origin_line_id)
            if base is None:
                new_pipes.append(l)
                continue
            # skip lines that were never sized
            if not hasattr(l, "inner_diameter_required"):
                continue
            if not hasattr(base, "inner_diameter_required"):
                continue
            
            if l.inner_diameter_required > base.inner_diameter_required:
                pipes_changed.append((base, l))
            else:
                l.dn = base.dn
                l.inner_diameter = base.inner_diameter
                l.outer_diameter = base.outer_diameter
                l.dn_index = base.dn_index
                l.length = min(l.length, base.length)

    scenario.dhn_pipe_changes = {
        "new_pipes": new_pipes,
        "pipes_changed": pipes_changed,
    }



def create_scenario ( 
        dhn_demand_path = None,
        dhn_nodes_path = None,
        elec_consumption_path =  None,
        elec_production_path =  None,
        elec_nodes_path = None,
        climate_file_path = None, 
        streets_path = None,
        random_seed: int = 42,
        mode: str = "DESIGN", 
        base = None,
        design = None
        ):

    if mode == "DESIGN":
        district_heating_network =  District_Heating_System(demand_path = dhn_demand_path,
                                         buildings_geo_path = dhn_nodes_path,
                                         streets_geo_path = streets_path,
                                         climate_file_path = climate_file_path,
                                         random_seed = random_seed)
        
        electrical_lv_network =      GridGenerator.from_geo(
                                        nodes_path= elec_nodes_path,
                                        streets_path= streets_path,
                                        building_consumption_csv= elec_consumption_path,
                                        building_production_csv= elec_production_path
                                    )
        
        district_heating_network.solve(time_frame = range(1,8760), #range in hour
                    return_pressure = 150_000.0, #Pa
                    min_pressure_difference = 100_000, #Pa
                    supply_temperature = 20, #C
                    required_temperature = 60, #C
                    ground_temperature_profile = 7.2, #C
                    solve_no_loss = False,
                    use_numba = True)
        
            
        for grid in electrical_lv_network.grids:
            grid.topology =analyze_topology(grid)
            grid.solve(timestep_hours=1.0)
            apply_results_to_grid(grid, grid.results)
        
        
        dhn_diagnoses = []
        # for _, system in district_heating_network.subsystems.items():
        #     dhn_diag = diagnose_dhn_faults(
        #         system,
        #         max_supply_heat=500_000, #W
        #         boiling_margin=0.10, 
        #         max_pipe_pressure=2_000_000, #Pa
        #         max_dp_per_m=200.0, #Pa/m
        #         max_pump_power=1_000, #W
        #     )
        #     dhn_diagnoses = dhn_diagnoses + dhn_diag.findings
        
            
        # gf = diag.generate_grid_fault_report_object(
        #     electrical_lv_network,
        #     min_voltage = 380,
        #     max_voltage = 420,
        #     max_supply_apparent_power = 1_000_000,
        # )
        elec_diagnoses = []
        # for i, fault in gf.grids.items():
        #     dg = diag.diagnose_grid_faults(fault,i)
        #     elec_diagnoses = elec_diagnoses + dg.findings
        
        
        
        Buildings = build_buildings_dict(electrical_lv_network, district_heating_network)
        Grids = electrical_lv_network.grids
        DistrictHeatingSystems = district_heating_network.subsystems
        
        S = Scenario(
            buildings = Buildings,
            Electrical_Network = Grids,
            District_Heating_Systems = DistrictHeatingSystems,
            Electrical_Diagnoses = elec_diagnoses,
            District_Heating_Diagnoses = dhn_diagnoses,
            )
        clean_scenario_fast(S)

        
        return S
    
    if mode == "BASELINE":

        scenario = copy.deepcopy(base)
        scenario = apply_dhn_changes(scenario = scenario, 
                                     csv_path = dhn_demand_path, 
                                     baseline_obj = scenario, 
                                     climate_path = climate_file_path)

        clean_scenario_fast(scenario)  
        return scenario
    
    
    if mode == "RETROFIT":

        scenario = copy.deepcopy(base)
        scenario = apply_dhn_changes(scenario = scenario, 
                                     csv_path = dhn_demand_path, 
                                     baseline_obj = scenario, 
                                     climate_path = climate_file_path)
        
        scenario = apply_grid_changes(
                            scenario,
                            consumption_csv_path= elec_consumption_path,
                            production_csv_path= elec_production_path,
                            baseline_obj = base)

        clean_scenario_fast(scenario)

            
        dhn_diagnoses = []
        for _, system in scenario.District_Heating_Systems.items():
            dhn_diag = diagnose_dhn_faults(
                system,
                max_supply_heat=500_000, #W
                boiling_margin=0.10, 
                max_pipe_pressure=2_000_000, #Pa
                max_dp_per_m=200.0, #Pa/m
                max_pump_power=1_000, #W
            )
            dhn_diagnoses = dhn_diagnoses + dhn_diag.findings
        scenario.District_Heating_Diagnoses = dhn_diagnoses
   
        gf = diag.generate_grid_fault_report_object(
            scenario.Electrical_Network,
            min_voltage = 380,
            max_voltage = 420,
            max_supply_apparent_power = 1_000_000,
        )
        elec_diagnoses = []
        for i, fault in gf.grids.items():
            dg = diag.diagnose_grid_faults(fault,i)
            elec_diagnoses = elec_diagnoses + dg.findings
        scenario.Electrical_Diagnoses = elec_diagnoses  

        return scenario

def build_buildings_dict(electrical_lv_network, district_heating_network):
    buildings = {}

    for node in electrical_lv_network.nodes.values():
        if getattr(node, "node_type", None) != "building":
            continue
        bid = node.building_id
        if bid not in buildings:
            buildings[bid] = Building(building_id=bid)
        buildings[bid].x = getattr(node, "x", None)
        buildings[bid].y = getattr(node, "y", None)
        consumption = getattr(node, "building_consumption", None)
        production = getattr(node, "building_production", None) 
        buildings[bid].electricity_bought = np.maximum(consumption - production, 0)
        buildings[bid].electricity_sold = np.maximum(production - consumption, 0)
        buildings[bid].voltage = getattr(node, "V", None)
    for dh in district_heating_network.subsystems.values():
        for node in dh.nodes:
            if getattr(node, "node_type", None) != "consumer":
                continue
            bid = getattr(node, "node_id", None)
            if bid is None:
                continue
            if bid not in buildings:
                buildings[bid] = Building(building_id=bid)
            buildings[bid].energy_bought_from_dhn = getattr(node, "demand", None)
            buildings[bid].dhn_supply_pressure = getattr(node, "supply_pressure", None)/100_000
            buildings[bid].dhn_supply_temperature = getattr(node, "supply_temperature", None)

    return buildings



def update_electrical_profiles(
    scenario,
    consumption_csv_path: str,
    production_csv_path: str,
):
    import pandas as pd
    import numpy as np
    cons_df = consumption_csv_path if isinstance(consumption_csv_path, pd.DataFrame) else pd.read_csv(consumption_csv_path, sep=";")
    prod_df = production_csv_path if isinstance(production_csv_path, pd.DataFrame) else pd.read_csv(production_csv_path, sep=";")



    cons_df.columns = cons_df.columns.astype(str)
    prod_df.columns = prod_df.columns.astype(str)

    for grid in scenario.Electrical_Network:
        for node in grid.nodes:
            if getattr(node, "node_type", None) != "building":
                continue

            bid = str(getattr(node, "building_id", None))
            if bid is None:
                continue

            if bid in cons_df.columns:
                node.building_consumption = cons_df[bid].to_numpy(dtype=float)
            else:
                node.building_consumption = np.zeros(
                    len(prod_df.index), dtype=float
                )

            if bid in prod_df.columns:
                node.building_production = prod_df[bid].to_numpy(dtype=float)
            else:
                node.building_production = np.zeros(
                    len(cons_df.index), dtype=float
                )
                
def recompute_grid_currents(
    scenario,
    *,
    slack_voltage: float = 400.0,
    timestep_hours: float = 1.0,
):
    from eureca_pubem.grid.GridTopology import analyze_topology
    from eureca_pubem.grid.powerflow_timeseries import apply_results_to_grid

    for grid in scenario.Electrical_Network:
        grid.topology = analyze_topology(grid)

        grid.solve(timestep_hours=timestep_hours)

        apply_results_to_grid(grid, grid.results)
        

def resize_grid_against_baseline(
    baseline_obj,
    scenario,
    *,
    cable_confidence_factor: float = 1.0,
):
    new_lines = []
    lines_changed = []

    baseline = baseline_obj
    for grid in baseline.Electrical_Network:
        ensure_origin_line_id_grid(grid)

    for grid in scenario.Electrical_Network:
        ensure_origin_line_id_grid(grid)

    baseline_lookup = {
        l.origin_line_id: l
        for grid in baseline.Electrical_Network
        for l in grid.lines
        if getattr(l, "origin_line_id", None) is not None
    }

    for grid in scenario.Electrical_Network:
        for l in grid.lines:
            oid = getattr(l, "origin_line_id", None)

            if oid is None:
                new_lines.append(l)
                continue

            base = baseline_lookup.get(oid)

            if base is None:
                new_lines.append(l)
                continue

            I_req = np.nanmax(getattr(l, "I", np.array([0.0])))
            I_base = np.nanmax(getattr(base, "I", np.array([0.0])))

            if I_req * cable_confidence_factor > getattr(base, "max_current_a", 0.0):
                lines_changed.append((base, l))
                upgrade_line_capacity(l, I_req, cable_confidence_factor)
            else:
                freeze_unchanged_lines(l, base)

    scenario.grid_line_changes = {
        "new_lines": new_lines,
        "lines_changed": lines_changed,
    }
    
    
def solve_electrical_network(
    scenario,
    *,
    slack_voltage: float = 400.0,
    timestep_hours: float = 1.0,
):
    from eureca_pubem.grid.GridTopology import analyze_topology
    from eureca_pubem.grid.powerflow_timeseries import apply_results_to_grid

    for grid in scenario.Electrical_Network:
        for node in grid.nodes:
            if getattr(node, "node_type", None) == "supply":
                node.set_voltage = slack_voltage

        grid.topology = analyze_topology(grid)

        grid.solve(timestep_hours=timestep_hours)

        apply_results_to_grid(grid, grid.results)

def ensure_origin_line_id_grid(grid):
    for l in grid.lines:
        if not hasattr(l, "origin_line_id") or l.origin_line_id is None:
            l.origin_line_id = l.line_id
            
            
def upgrade_line_capacity(line, I_required, cable_confidence_factor=1.0):
    import json
    from pathlib import Path
    import numpy as np

    base_path = Path(__file__).resolve().parent.parent / "grid/standard_cables_lv.json"
    with open(base_path, "r") as f:
        catalog = json.load(f)

    I_target = float(I_required) * float(cable_confidence_factor)

    candidates = sorted(
        catalog,
        key=lambda x: float(x.get("max_current_a", 0.0))
    )

    for c in candidates:
        if float(c.get("max_current_a", 0.0)) >= I_target:
            line.cable_type = int(c["id"])
            line.r_ohm_per_km = float(c["r_ohm_per_km"])
            line.x_ohm_per_km = float(c["x_ohm_per_km"])
            line.max_current_a = float(c["max_current_a"])
            return

    c = candidates[-1]
    line.cable_type = int(c["id"])
    line.r_ohm_per_km = float(c["r_ohm_per_km"])
    line.x_ohm_per_km = float(c["x_ohm_per_km"])
    line.max_current_a = float(c["max_current_a"])
    
def freeze_unchanged_lines(line, baseline_line):
    attrs = [
        "cable_type",
        "r_ohm_per_km",
        "x_ohm_per_km",
        "max_current_a",
    ]

    for a in attrs:
        if hasattr(baseline_line, a):
            setattr(line, a, getattr(baseline_line, a))    
    