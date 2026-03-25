from __future__ import annotations

from dataclasses import dataclass 
from typing import Dict, List, Set
import numpy as np 
import networkx as nx
import pandas as pd 
import geopandas as gpd
import networkx as nxbuild_dhn_from_geometry
import eureca_pubem.dhn.DH_tools as toolbox
import eureca_pubem.dhn.DH_tools_fast as fast_toolbox
from eureca_pubem.dhn.DH_design import build_dhn_from_geometry, parse_buildings_nodes
from eureca_pubem.dhn.DH_objects import Node, Line, DistrictHeating
from eureca_pubem.dhn.DH_plots import plot_dhn_capacity_layout
from time import time

class District_Heating_System:
    REQUIRED_SHEETS = ("demands", "supply", "lines")
    REQUIRED_SHEETS_GEO = ("demands",)

    def __init__(self,
                 demand_path: str,
                 random_seed,
                 buildings_geo_path: str = None,
                 streets_geo_path: str = None,
                 climate_file_path: str = None):
    

    
        use_geo = buildings_geo_path is not None and streets_geo_path is not None
        if use_geo:
            self.filepath = demand_path

            demands_df = demand_path if isinstance(demand_path, pd.DataFrame) else pd.read_csv(demand_path, sep=";")
            demand_node_ids = self._parse_demand_headers(demands_df)
            demand_matrix = demands_df.to_numpy(dtype=float)
            n_timesteps = demand_matrix.shape[0]
            self.street_path = streets_geo_path
            nodes, _ = parse_buildings_nodes(buildings_geo_path, streets_geo_path)
            self.prenodes = nodes
            supplies = [n for n in nodes.values() if n.node_type == "supply"]
            has_capacity = all(hasattr(n, "capacity") for n in supplies)
            if any(hasattr(n, "capacity") for n in supplies) and not has_capacity:
                raise ValueError("Either all supply nodes must define capacity, or none.")
                
                
            if has_capacity:
                
                for nid, n in nodes.items():
                    if n.node_type == "consumer":
                        if nid not in demand_node_ids:
                            raise ValueError(f"No demand column for consumer {nid}")
                        j = demand_node_ids.index(nid)
                        n.peak_demand = float(np.max(demand_matrix[:, j]))
                assignments = toolbox.assign_consumers_to_supplies(
                    nodes,
                    concurrency_fn=lambda n: (65.0 * (1.2 ** (-0.05 * n)) + 35.0) / 100.0,
                    fail_on_unassigned=False,
                )

                self.assignments = assignments
                def load_gdf(obj):
                    if isinstance(obj, (str, bytes)):
                        return gpd.read_file(obj)
                    elif isinstance(obj, gpd.GeoDataFrame):
                        return obj
                    else:
                        raise TypeError("Input must be a GeoDataFrame or a path to a file.")    
                buildings_gdf = load_gdf(buildings_geo_path)
                
                self.subsystems = {}
                
                for supply_id, consumer_ids in assignments.items():
                    keep_ids = set([supply_id] + consumer_ids)
                
                    buildings_subset = buildings_gdf[
                        buildings_gdf["id"].astype(int).isin(keep_ids)
                    ]
                
                    if buildings_subset.empty:
                        continue
                    # self.dhn = build_dhn_from_geometry(
                    #     demand_csv_path = demand_path,
                    #     streets_path=streets_geo_path,
                    #     buildings_path=buildings_subset
                    # )

                    geo_nodes, geo_lines = build_dhn_from_geometry(
                        demand_csv_path = demand_path,
                        streets_path=streets_geo_path,
                        buildings_path=buildings_subset
                    )

                    
                    nodes_map = {}
                    for n in geo_nodes:
                        if n.node_id in demand_node_ids:
                            j = demand_node_ids.index(n.node_id)
                            n.demand = demand_matrix[:, j]
                        else:
                            n.demand = np.zeros(n_timesteps)
                
                        if n.node_type == "consumer":
                            n.peak_demand = float(np.max(n.demand))
                
                        nodes_map[n.node_id] = n
                
                    lines_map = {l.line_id: l for l in geo_lines}
                
                    name = f"District Heating Supply {supply_id}"
                    system = DistrictHeating(
                        name,
                        list(nodes_map.values()),
                        list(lines_map.values()),
                    )
                    system.node_ids = list(nodes_map.keys())
                    self.subsystems[name] = system
            else:

                geo_nodes, geo_lines = build_dhn_from_geometry(
                    streets_path=streets_geo_path,
                    buildings_path=buildings_geo_path,
                )
        
                self.nodes = {}
                for n in geo_nodes:
                    if n.node_id in demand_node_ids:
                        j = demand_node_ids.index(n.node_id)
                        n.demand = demand_matrix[:, j]
                    else:
                        n.demand = np.zeros(n_timesteps)
        
                    self.nodes[n.node_id] = n
        
                self.lines = {l.line_id: l for l in geo_lines}
        
                all_nodes = set(self.nodes.keys())
                components = self._connected_components(all_nodes, geo_lines)
                for n in self.nodes.values():
                    if n.node_type == "consumer":
                        n.peak_demand = float(np.max(n.demand))
                
            

              

        else:
            self.filepath = demand_path

            xls = pd.ExcelFile(demand_path)
        
            demands_df = pd.read_excel(demand_path, sheet_name="demands", header=0)
            demand_node_ids = self._parse_demand_headers(demands_df)
            demand_matrix = demands_df.to_numpy(dtype=float)
            n_timesteps = demand_matrix.shape[0]
            self.street_path = streets_geo_path
            missing = [s for s in self.REQUIRED_SHEETS if s not in xls.sheet_names]
            if missing:
                raise ValueError(f"Missing required sheet(s): {missing}")
    
            supply_df = pd.read_excel(demand_path, sheet_name="supply", header=0)
            lines_df = pd.read_excel(demand_path, sheet_name="lines", header=0)
    
            supply_nodes = self._parse_supply_nodes(supply_df)
            lines = self._parse_lines(lines_df)
    
            all_line_nodes = {l.start_node for l in lines} | {l.end_node for l in lines}
            all_nodes = set(demand_node_ids) | all_line_nodes
    
            self._validate_supply_nodes(supply_nodes, all_nodes)
            self._validate_lines(lines, all_nodes)
    
            self.nodes = {}
            for nid in sorted(all_nodes):
                if nid in demand_node_ids:
                    j = demand_node_ids.index(nid)
                    demand = demand_matrix[:, j]
                else:
                    demand = np.zeros(n_timesteps)
    
                if nid in supply_nodes:
                    ntype = "supply"
                elif np.any(demand != 0):
                    ntype = "consumer"
                else:
                    ntype = "connection"
    
                node = Node(nid, ntype)
                node.demand = demand
                self.nodes[nid] = node
    
            self.lines = {l.line_id: l for l in lines}
            components = self._connected_components(all_nodes, lines)
    
        if has_capacity == False:
            self.subsystems = {}
            for i, comp in enumerate(components, start=1):
                name = f"District Heating {i}"
                comp_nodes = [self.nodes[n] for n in sorted(comp)]
                comp_lines = [
                    l for l in self.lines.values()
                    if l.start_node in comp and l.end_node in comp
                ]
                self.subsystems[name] = DistrictHeating(name, comp_nodes, comp_lines)
                system = self.subsystems[name]
                system.node_ids = [n.node_id for n in system.nodes]

        self.create_matrices(random_seed)
    
        for _, system in self.subsystems.items():
            toolbox.lines_sizing(system)
    
        self.precompute_cache(
            epw_file=climate_file_path,
            depth=4,
            minimum_flow_velocity=0.1,
            required_node_supply_temperature=60,
            expected_node_return_temperature=40,
            seed=random_seed,
        )


    def _parse_demand_headers(self, df):
        node_ids = []
        for c in df.columns:
            try:
                node_ids.append(int(c))
            except ValueError:
                continue
        if len(set(node_ids)) != len(node_ids):
            raise ValueError("Duplicate node ids in demands sheet.")
        return node_ids

    def _parse_supply_nodes(self, df: pd.DataFrame) -> Set[int]:
        if "supply_nodes" not in df.columns:
            raise ValueError("Sheet 'supply' must contain column 'supply_nodes'.")
        return set(int(v) for v in df["supply_nodes"].dropna())

    def _parse_lines(self, df: pd.DataFrame) -> List[Line]:
        required = {
            "line_id",
            "start_node",
            "end_node",
            "length",
            "inner_diameter",
            "outer_diameter",
            "pipe_conductivity",
            "insulation_thickness",
            "insulation_conductivity",
            "roughness",
            "pipe_specific_heat_capacity",
        }
        if not required.issubset(df.columns):
            raise ValueError(f"Sheet 'lines' must contain columns {required}")

        lines = []
        for _, r in df.iterrows():
            lines.append(
                Line(
                    int(r["line_id"]),
                    int(r["start_node"]),
                    int(r["end_node"]),
                    float(r["length"]),
                    float(r["inner_diameter"]),
                    float(r["outer_diameter"]),
                    float(r["pipe_conductivity"]),
                    float(r["insulation_conductivity"]),
                    float(r["insulation_thickness"]),
                    float(r["roughness"]),
                    float(r["pipe_specific_heat_capacity"]),
                    0.0,
                    0.0,
                    float(r["max_pressure"]) if "max_pressure" in df.columns and not pd.isna(r["max_pressure"]) else None
                )
            )
        if len({l.line_id for l in lines}) != len(lines):
            raise ValueError("Duplicate line_id in lines sheet.")
        return lines

    def _validate_supply_nodes(self, supply: Set[int], nodes: Set[int]):
        missing = supply - nodes
        if missing:
            raise ValueError(f"Supply nodes not found in network: {missing}")

    def _validate_lines(self, lines: List[Line], nodes: Set[int]):
        for l in lines:
            if l.start_node not in nodes or l.end_node not in nodes:
                raise ValueError(f"Line {l.line_id} references unknown node.")

    def _connected_components(self, node_ids: Set[int], lines: List[Line]) -> List[Set[int]]:
        adj: Dict[int, Set[int]] = {n: set() for n in node_ids}
        for l in lines:
            adj[l.start_node].add(l.end_node)
            adj[l.end_node].add(l.start_node)

        visited = set()
        components = []

        for n in node_ids:
            if n in visited:
                continue
            stack = [n]
            comp = set()
            while stack:
                x = stack.pop()
                if x in visited:
                    continue
                visited.add(x)
                comp.add(x)
                stack.extend(adj[x] - visited)
            components.append(comp)

        return components
    
    def create_matrices (self , seed = 17):
        for _, system in self.subsystems.items():
            toolbox.incidence_matrix(system)

            toolbox.cutset_matrix_with_supply_datum(system, seed=seed)
            toolbox.loop_matrix(system)
            
    def precompute_cache (self, 
                          epw_file,
                          depth,
                          minimum_flow_velocity = 0.15,
                          required_node_supply_temperature = 60, 
                          expected_node_return_temperature = 40,
                          seed = 17):
        for _, system in self.subsystems.items():
            toolbox.initialize_mass_flow_guesses_tree_routing_timeseries(system,
                                                                         min_velocity= minimum_flow_velocity,
                                                                         Tin_default=required_node_supply_temperature, 
                                                                         Tout_default=expected_node_return_temperature)
            toolbox.build_tree_pressure_cache(system, seed=seed)
            toolbox.precompute_line_UA(system)
            system.ground_temperature_array = toolbox.earth_temperature_from_epw(epw_file,
                                                                               depth)
            
    def solve (self,
               time_frame = range(1,8760),
               return_pressure = 150_000.0,
               min_pressure_difference = 100_000,
               supply_temperature = 80,
               required_temperature = 60,
               ground_temperature_profile = 7.2,
               solve_no_loss = False,
               use_numba = False):

            
        for _, system in self.subsystems.items():

            if use_numba:

                fast_toolbox.solve_timeseries_fast_tree(system, time_frame, 
                                                        return_pressure_default= return_pressure,
                                                        required_pressure_difference_default=min_pressure_difference,
                                                        supply_temperature=supply_temperature,
                                                        min_supply_temperature=required_temperature)
                if solve_no_loss:
                    for t in time_frame:
                        toolbox.solve_thermal_matrix_meanT_tree_no_loss(system, t)
            else:
                for t in time_frame:
                    toolbox.compute_pressures_fast_tree(system, t)
                    toolbox.solve_thermal_matrix_meanT_tree(system, t)
            if solve_no_loss:
                toolbox.solve_thermal_matrix_meanT_tree_no_loss(system, t)
                
            toolbox.add_line_thermal_losses(system)
            toolbox.add_hourly_energy_and_pump_metrics(system)

    def plot(self):
        for name, system in self.subsystems.items():
            toolbox.plot_info(system, name)
    def plot_layout_subs(self)   :     
        plot_dhn_capacity_layout(
            self.subsystems,
            streets_gdf=gpd.read_file(self.street_path),
            unassigned_nodes=[
                n for n in self.prenodes.values()
                if n.node_type == "consumer" and n.assigned_supply_id is None
            ],
            save_path="dhn_capacity_layout.svg",
        )
        
    @classmethod
    def from_subsystems(cls,
                        subsystems: Dict[str, DistrictHeating],
                        random_seed: int,
                        climate_file_path: str = None):
        
        obj = cls.__new__(cls)  # bypass __init__
        
        obj.subsystems = subsystems

        obj.filepath = None
        obj.street_path = None
        for _, system in obj.subsystems.items():
            obj._clean_subsystem_topology(system)



        obj.create_matrices(random_seed)
    
        # sizing
        for _, system in obj.subsystems.items():
            toolbox.lines_sizing(system)
    
        # precompute
        obj.precompute_cache(
            epw_file=climate_file_path,
            depth=4,
            minimum_flow_velocity=0.1,
            required_node_supply_temperature=60,
            expected_node_return_temperature=40,
            seed=random_seed,
        )
    
        return obj


    
    def _clean_subsystem_topology(self, system: DistrictHeating):
        """
        Keep only active nodes/lines and enforce a single connected tree.
        Remove isolated nodes.
        """
    
        # --- Filter active nodes
        active_nodes = {
            n.node_id: n
            for n in system.nodes
            if getattr(n, "active", True)
        }
    
        # --- Filter active lines with valid endpoints
        active_lines = [
            l for l in system.lines
            if getattr(l, "active", True)
            and l.start_node in active_nodes
            and l.end_node in active_nodes
        ]
    
        # --- Build graph
        G = nx.Graph()
    
        for l in active_lines:
            G.add_edge(l.start_node, l.end_node, object=l)
    
        if G.number_of_nodes() == 0:
            system.nodes = []
            system.lines = []
            system.node_ids = []
            return
    
        # --- Extract spanning tree (forest-safe)
        T = nx.minimum_spanning_tree(G)
    
        kept_lines = []
        connected_nodes = set()
    
        for u, v in T.edges:
            edge_data = G.get_edge_data(u, v)
            kept_lines.append(edge_data["object"])
            connected_nodes.add(u)
            connected_nodes.add(v)
    
        # --- Keep only nodes actually in tree
        system.nodes = [
            active_nodes[nid]
            for nid in connected_nodes
        ]
    
        system.lines = kept_lines
        system.node_ids = list(connected_nodes)

