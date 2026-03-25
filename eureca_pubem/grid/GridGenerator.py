from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
import math 
from eureca_pubem.grid.Grid_design import build_grid_from_geometry, populate_grid_system
from eureca_pubem.grid.Grid_classes import GridSystem, Grid, Node, Line, Inverter


class GridGenerator:

    @classmethod
    def from_geo(
        cls,
        *,
        nodes_path,
        streets_path,
        building_consumption_csv,
        building_production_csv,
        pf_load: float = 0.95,
        load_lagging: bool = True,
        pf_pv: float = 1.0,
        pv_lagging: bool = False,
        voltage_lv_min: float =380.0,
        cable_confidence_factor: float = 1.0,
    ) -> GridSystem:
        
        
        GS =  build_grid_from_geometry(
        nodes_path=nodes_path,
        streets_path=streets_path,
        )
        GS_populated =  populate_grid_system(
            grid = GS,
            consumption_path= building_consumption_csv,
            production_path= building_production_csv,
            voltage_lv= voltage_lv_min,
            confidence_factor = cable_confidence_factor
        )
        return GS_populated

    @classmethod
    def from_excel(
        cls,
        *,
        excel_path: str,
        nodes_sheet: str = "buses",
        lines_sheet: str = "lines",
        linetypes_sheet: str = "linetypes",
        substation_sheet: str = "substation",
        building_consumption_sheet: str = "building_consumption",
        building_production_sheet: str = "building_production",
        pf_load: float = 0.95,
        load_lagging: bool = True,
        pf_pv: float = 1.0,
        pv_lagging: bool = False,
    ) -> GridSystem:
        nodes_df = pd.read_excel(excel_path, sheet_name=nodes_sheet)
        lines_df = pd.read_excel(excel_path, sheet_name=lines_sheet)
        linetypes_df = pd.read_excel(excel_path, sheet_name=linetypes_sheet)
        substation_df = pd.read_excel(excel_path, sheet_name=substation_sheet)
    
        cons_df, n_cons = cls._read_timeseries_sheet(excel_path, building_consumption_sheet)
        prod_df, n_prod = cls._read_timeseries_sheet(excel_path, building_production_sheet)
    
        if n_cons != n_prod:
            raise ValueError("Consumption and production length mismatch")
    
        cable_lookup = cls._build_linetype_lookup(linetypes_df)
        sub_lookup = cls._build_substation_lookup(substation_df)
    
        nodes: Dict[int, Node] = {}
    
        for _, r in nodes_df.iterrows():
            bus_id = int(r["Bus ID"])
            node_type = str(r["Type"]).strip().lower()
            x = float(r["X"])
            y = float(r["Y"])
    
            building_id = None if pd.isna(r.get("Building ID")) else int(r["Building ID"])
    
            max_power = set_voltage = set_voltage_angle = None
            if node_type == "supply":
                max_power, set_voltage, set_voltage_angle = sub_lookup[bus_id]
    
            bcons = bprod = None
            if node_type == "building":
                col = str(building_id)
                bcons = cons_df[col].to_numpy(dtype=float, copy=True) if col in cons_df.columns else np.zeros(n_cons)
                bprod = prod_df[col].to_numpy(dtype=float, copy=True) if col in prod_df.columns else np.zeros(n_cons)
    
            nodes[bus_id] = Node(
                bus_id=bus_id,
                node_type=node_type,
                x=x,
                y=y,
                building_id=building_id,
                max_power=max_power,
                set_voltage=set_voltage,
                set_voltage_angle=set_voltage_angle,
                building_consumption=bcons,
                building_production=bprod,
                inverter=None,
                pf_load=pf_load,
                load_lagging=load_lagging,
                pf_pv=pf_pv,
                pv_lagging=pv_lagging,
            )
    
        lines: List[Line] = []
    
        for _, r in lines_df.iterrows():
            line_id = int(r["Line ID"])
            from_bus = int(r["From bus"])
            to_bus = int(r["To bus"])
            cable_type = None if pd.isna(r["Cable type"]) else int(r["Cable type"])
    
            n1 = nodes[from_bus]
            n2 = nodes[to_bus]
            length_m = math.hypot(n1.x - n2.x, n1.y - n2.y)
    
            r_ohm = x_ohm = max_i = None
            if cable_type is not None:
                r_ohm, x_ohm, max_i = cable_lookup[cable_type]
    
            lines.append(
                Line(
                    line_id=line_id,
                    from_bus=from_bus,
                    to_bus=to_bus,
                    cable_type=cable_type,
                    length_m=length_m,
                    r_ohm_per_km=r_ohm,
                    x_ohm_per_km=x_ohm,
                    max_current_a=max_i,
                )
            )
    
        return GridSystem(nodes=nodes, lines=lines)

    @staticmethod
    def _read_timeseries_sheet(
        excel_path: str,
        sheet_name: str,
    ) -> Tuple[pd.DataFrame, int]:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
        df.columns = [str(c).strip() for c in df.columns]
        return df, len(df)

    @staticmethod
    def _build_linetype_lookup(
        linetypes_df: pd.DataFrame,
    ) -> Dict[int, Tuple[float, float, float]]:
        out = {}
        for _, r in linetypes_df.iterrows():
            cid = int(r["Cable ID"])
            out[cid] = (
                float(r["Resistance"]),
                float(r["Reactance"]),
                float(r["Max current"]),
            )
        return out

    @staticmethod
    def _build_substation_lookup(
        substation_df: pd.DataFrame,
    ) -> Dict[int, Tuple[float, float, float]]:
        out = {}
        for _, r in substation_df.iterrows():
            bid = int(r["bus id"])
            out[bid] = (
                float(r["Max Power"]),
                float(r["Set_voltage"]),
                float(r["Set_voltage_angle"]),
            )
        return out

    def attach_inverters_for_producers(
        self,
        *,
        capacity_factor: float = 1.2,
        default_efficiency: float = 0.98,
        overwrite_existing: bool = False,
    ) -> None:
        for bus_id, node in list(self.nodes.items()):
            if node.node_type != "building":
                continue
    
            prod = node.building_production
            if prod is None or prod.size == 0:
                continue
    
            max_prod = float(np.nanmax(prod))
            if max_prod <= 0:
                continue
    
            if node.inverter is not None and not overwrite_existing:
                continue
    
            inv = Inverter(
                capacity_w=capacity_factor * max_prod,
                efficiency=default_efficiency,
            )
    
            self.nodes[bus_id] = Node(
                bus_id=node.bus_id,
                node_type=node.node_type,
                x=node.x,
                y=node.y,
                building_id=node.building_id,
                max_power=node.max_power,
                set_voltage=node.set_voltage,
                set_voltage_angle=node.set_voltage_angle,
                building_consumption=node.building_consumption,
                building_production=node.building_production,
                inverter=inv,
                pf_load=node.pf_load,
                load_lagging=node.load_lagging,
                pf_pv=node.pf_pv,
                pv_lagging=node.pv_lagging,
            )
