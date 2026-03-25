import geopandas as gpd
import numpy as np
import networkx as nx
import pandas as pd
from typing import  Optional
import json
from pathlib import Path

from shapely.geometry import  LineString, MultiLineString
from shapely.ops import split

from eureca_pubem.grid.Grid_classes import Node, Line, GridSystem, Grid

def snap_components(graph):

    components = list(nx.connected_components(graph))
    main = max(components, key=len)
    main_nodes = list(main)

    for comp in components:
        if comp == main:
            continue

        comp_nodes = list(comp)

        best_u = None
        best_v = None
        best_d = float("inf")

        for u in comp_nodes:
            for v in main_nodes:
                d = math.dist(u, v)
                if d < best_d:
                    best_d = d
                    best_u = u
                    best_v = v

        graph.add_edge(best_u, best_v, weight=best_d, kind="snap")

    return graph

def build_grid_from_geometry(
    nodes_path: str,
    streets_path: str,
    *,
    crs_metric = "EPSG:32633",
    grid_resolution=None,
    grid_penalty=None,
) -> GridSystem:
    def load_gdf(obj):
        if isinstance(obj, (str, bytes)):
            return gpd.read_file(obj)
        elif isinstance(obj, gpd.GeoDataFrame):
            return obj
        else:
            raise TypeError("Input must be a GeoDataFrame or a path to a file.")    
    nodes_gdf = load_gdf(nodes_path)
    def load_gdf(obj):
        if isinstance(obj, (str, bytes)):
            return gpd.read_file(obj)
        elif isinstance(obj, gpd.GeoDataFrame):
            return obj
        else:
            raise TypeError("Input must be a GeoDataFrame or a path to a file.") 
    streets_gdf = load_gdf(streets_path)

    if nodes_gdf.crs is None:
        raise ValueError("nodes_gdf has no CRS")
    if streets_gdf.crs is None:
        raise ValueError("streets_gdf has no CRS")
    if crs_metric is None:
        raise ValueError("crs_metric must be provided")

    nodes_gdf = nodes_gdf.to_crs(crs_metric)
    streets_gdf = streets_gdf.to_crs(crs_metric)

    if nodes_gdf.crs != streets_gdf.crs:
        raise ValueError("nodes and streets CRS mismatch after reprojection")

    streets = sanitize_street_geometry(streets_gdf, crs_metric)
    street_graph = build_street_graph(streets)
    grid_graph = build_auxiliary_grid(streets, grid_resolution, grid_penalty)
    routing_graph = merge_graphs(street_graph, grid_graph)
    routing_graph = attach_terminal_nodes(routing_graph, nodes_gdf)
    assignment = assign_buildings_to_supplies(routing_graph, nodes_gdf)
    routing_graph = snap_components(routing_graph)
    trees = extract_supply_trees(routing_graph, nodes_gdf, assignment)

    return realize_grid_system(trees, nodes_gdf)




def sanitize_street_geometry(
    streets_gdf: gpd.GeoDataFrame,
    crs_metric = "EPSG:32633",
    *,
    tol: float = 1e-6,
) -> gpd.GeoDataFrame:
    if crs_metric is None:
        raise ValueError("crs_metric must be provided and must be a projected CRS")

    if streets_gdf.crs is None:
        raise ValueError("streets_gdf has no CRS")

    gdf = streets_gdf.to_crs(crs_metric)
    crs = gdf.crs

    gdf["geometry"] = gdf.geometry.apply(
        lambda g: list(g.geoms) if isinstance(g, MultiLineString) else [g]
    )
    gdf = gdf.explode("geometry", ignore_index=True)

    gdf = gdf[gdf.geometry.apply(lambda g: g.length) > tol]

    nd = int(abs(np.log10(tol)))
    gdf["geometry"] = gdf.geometry.apply(
        lambda g: LineString(np.round(np.asarray(g.coords), nd))
    )

    gdf = gdf.drop_duplicates(subset="geometry").reset_index(drop=True)

    lines = list(gdf.geometry)
    split_lines = []

    for i, li in enumerate(lines):
        geom = li
        for j, lj in enumerate(lines):
            if j <= i:
                continue
            if not geom.intersects(lj):
                continue
            ip = geom.intersection(lj)
            if ip.is_empty:
                continue
            try:
                geom = split(geom, ip)
            except Exception:
                continue
        if isinstance(geom, MultiLineString):
            split_lines.extend(list(geom.geoms))
        else:
            split_lines.append(geom)

    out = gpd.GeoDataFrame(geometry=split_lines, crs=crs)
    out = out[out.geometry.apply(lambda g: g.length) > tol].reset_index(drop=True)

    return out





def build_street_graph(
    streets_gdf: gpd.GeoDataFrame,
) -> nx.Graph:
    G = nx.Graph()

    for geom in streets_gdf.geometry:
        if geom is None or geom.is_empty:
            continue

        stack = [geom]

        while stack:
            g = stack.pop()

            if isinstance(g, MultiLineString):
                stack.extend(list(g.geoms))
                continue

            if not isinstance(g, LineString):
                continue

            coords = list(g.coords)
            if len(coords) < 2:
                continue

            for i in range(len(coords) - 1):
                u = (float(coords[i][0]), float(coords[i][1]))
                v = (float(coords[i + 1][0]), float(coords[i + 1][1]))

                if u not in G:
                    G.add_node(u)
                if v not in G:
                    G.add_node(v)

                w = np.hypot(v[0] - u[0], v[1] - u[1])

                if G.has_edge(u, v):
                    if w < G[u][v]["weight"]:
                        G[u][v]["weight"] = w
                else:
                    G.add_edge(u, v, weight=w, kind="street")

    return G




def build_auxiliary_grid(
    streets_gdf: gpd.GeoDataFrame,
    grid_resolution: Optional[float],
    grid_penalty: Optional[float],
) -> nx.Graph:
    bounds = streets_gdf.total_bounds
    xmin, ymin, xmax, ymax = bounds

    if grid_resolution is None:
        grid_resolution = max(xmax - xmin, ymax - ymin) / 30.0

    if grid_penalty is None:
        grid_penalty = 20.0

    xs = np.arange(xmin, xmax + grid_resolution, grid_resolution)
    ys = np.arange(ymin, ymax + grid_resolution, grid_resolution)

    G = nx.Graph()

    for x in xs:
        for y in ys:
            G.add_node((float(x), float(y)))

    for x in xs:
        for y in ys:
            u = (float(x), float(y))
            for dx, dy in (
                (grid_resolution, 0.0),
                (0.0, grid_resolution),
                (grid_resolution, grid_resolution),
                (grid_resolution, -grid_resolution),
            ):
                v = (x + dx, y + dy)
                if v[0] < xmin or v[0] > xmax or v[1] < ymin or v[1] > ymax:
                    continue
                v = (float(v[0]), float(v[1]))
                if v not in G:
                    continue
                w = grid_penalty * np.hypot(dx, dy)
                G.add_edge(u, v, weight=w, kind="grid")

    return G




def _nearest_node(point, nodes):
    d = nodes - point
    dist = np.hypot(d[:, 0], d[:, 1])
    i = int(np.argmin(dist))
    return tuple(nodes[i]), float(dist[i])


def merge_graphs(
    street_graph: nx.Graph,
    grid_graph: nx.Graph,
) -> nx.Graph:
    G = nx.Graph()

    for u, data in street_graph.nodes(data=True):
        G.add_node(u, **data)

    for u, v, data in street_graph.edges(data=True):
        G.add_edge(u, v, **data)

    for u, data in grid_graph.nodes(data=True):
        if u not in G:
            G.add_node(u, **data)

    for u, v, data in grid_graph.edges(data=True):
        if G.has_edge(u, v):
            if data.get("weight", np.inf) < G[u][v].get("weight", np.inf):
                G[u][v].update(data)
        else:
            G.add_edge(u, v, **data)

    street_nodes = np.array(list(street_graph.nodes), dtype=float)
    grid_nodes = np.array(list(grid_graph.nodes), dtype=float)

    for sn in street_graph.nodes:
        gn, d = _nearest_node(np.array(sn), grid_nodes)
        if not G.has_edge(sn, gn):
            G.add_edge(sn, gn, weight=d, kind="street_grid_link")

    return G



def attach_terminal_nodes(
    graph: nx.Graph,
    nodes_gdf: gpd.GeoDataFrame,
) -> nx.Graph:
    G = graph.copy()

    graph_nodes = np.array(list(G.nodes), dtype=float)
    if graph_nodes.size == 0:
        raise ValueError("graph has no nodes")

    for _, row in nodes_gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        p = np.array([float(geom.x), float(geom.y)])

        d = graph_nodes - p
        dist = np.hypot(d[:, 0], d[:, 1])
        idx = int(np.argmin(dist))

        nearest = tuple(graph_nodes[idx])
        u = (float(p[0]), float(p[1]))

        if u not in G:
            G.add_node(u)

        w = float(dist[idx])
        G.add_edge(u, nearest, weight=w, kind="terminal")

    return G


def assign_buildings_to_supplies(
    graph: nx.Graph,
    nodes_gdf: gpd.GeoDataFrame,
) -> dict[int, int]:
    supplies = {}
    buildings = {}

    for _, row in nodes_gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        p = np.array([float(geom.x), float(geom.y)])
        node_id = int(row["id"])
        if int(row["type"]) == 1:
            supplies[node_id] = p
        elif int(row["type"]) == -1:
            buildings[node_id] = p

    if not supplies:
        raise ValueError("no supply nodes found")

    assignment = {}

    for bid, bp in buildings.items():
        dmin = np.inf
        sid_min = None
        for sid, sp in supplies.items():
            d = np.hypot(bp[0] - sp[0], bp[1] - sp[1])
            if d < dmin:
                dmin = d
                sid_min = sid
        if sid_min is None:
            raise ValueError(f"building {bid} not assigned")
        assignment[bid] = sid_min

    return assignment




import networkx as nx
import geopandas as gpd
import math

def extract_supply_trees(graph: nx.Graph, nodes_gdf: gpd.GeoDataFrame, assignment: dict[int, int]) -> dict[int, nx.Graph]:
    nodes = list(graph.nodes)

    


    def snap(p):
        if p in graph:
            return p
        best = None
        best_d = float("inf")
        x, y = p
        for n in nodes:
            d = (n[0]-x)**2 + (n[1]-y)**2
            if d < best_d:
                best_d = d
                best = n
        return best

    geom_lookup = nodes_gdf.set_index("id").geometry

    supplies = {}
    buildings_by_supply = {}

    for _, row in nodes_gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        p = snap((float(geom.x), float(geom.y)))
        if int(row["type"]) == 1:
            supplies[int(row["id"])] = p

    for bid, sid in assignment.items():
        geom = geom_lookup[bid]
        p = snap((float(geom.x), float(geom.y)))
        buildings_by_supply.setdefault(sid, []).append(p)

    trees = {}

    for sid, sp in supplies.items():
        Gs = nx.Graph()
        targets = buildings_by_supply.get(sid, [])

        for tp in targets:
            try:
                path = nx.shortest_path(graph, source=sp, target=tp, weight="weight")
            except nx.NetworkXNoPath:
                continue

            for u, v in zip(path[:-1], path[1:]):
                data = graph[u][v]

                if not Gs.has_node(u):
                    Gs.add_node(u)
                if not Gs.has_node(v):
                    Gs.add_node(v)

                if not Gs.has_edge(u, v):
                    Gs.add_edge(u, v, weight=data["weight"], kind=data.get("kind"))

        trees[sid] = nx.minimum_spanning_tree(Gs, weight="weight")

    return trees


def realize_grid_system(
    trees: dict[int, nx.Graph],
    nodes_gdf: gpd.GeoDataFrame,
) -> GridSystem:
    nodes: dict[int, Node] = {}
    lines: list[Line] = []

    coord_to_bus: dict[tuple[float, float], int] = {}
    next_bus_id = 1
    next_line_id = 1

    def get_bus_id(coord):
        nonlocal next_bus_id
        if coord not in coord_to_bus:
            coord_to_bus[coord] = next_bus_id
            next_bus_id += 1
        return coord_to_bus[coord]

    terminal_map: dict[tuple[float, float], dict] = {}
    for _, row in nodes_gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        p = (float(geom.x), float(geom.y))
        terminal_map[p] = row

    for sid, G in trees.items():
        for coord in G.nodes:
            bus_id = get_bus_id(coord)

            if bus_id in nodes:
                continue

            row = terminal_map.get(coord, None)
            if row is not None:
                if int(row["type"]) == 1:
                    n = Node(
                        bus_id=bus_id,
                        x=coord[0],
                        y=coord[1],
                        node_type="supply",
                    )
                else:
                    n = Node(
                        bus_id=bus_id,
                        x=coord[0],
                        y=coord[1],
                        node_type="building",
                        building_id=int(row["id"]),
                    )
            else:
                n = Node(
                    bus_id=bus_id,
                    x=coord[0],
                    y=coord[1],
                    node_type="connection",
                )

            nodes[bus_id] = n

        for u, v, data in G.edges(data=True):
            bu = get_bus_id(u)
            bv = get_bus_id(v)

            length = float(np.hypot(v[0] - u[0], v[1] - u[1]))

            l = Line(
                line_id=next_line_id,
                from_bus=bu,
                to_bus=bv,
                length_m=length,
            )
            next_line_id += 1
            lines.append(l)
            
    return split_grids_from_trees(GridSystem(nodes=nodes, lines=lines), trees)

def split_grids_from_trees(
    grid_system: GridSystem,
    trees: dict[int, nx.DiGraph],
) -> list[Grid]:

    coord_to_bus = {
        (float(n.x), float(n.y)): n.bus_id
        for n in grid_system.nodes.values()
    }

    grids = []

    for _, T in trees.items():
        bus_ids = {coord_to_bus[c] for c in T.nodes if c in coord_to_bus}

        nodes = [grid_system.nodes[bid] for bid in bus_ids]

        lines = []
        for u, v in T.edges:
            bu = coord_to_bus.get(u)
            bv = coord_to_bus.get(v)
            if bu is None or bv is None:
                continue

            for l in grid_system.lines:
                if (
                    (l.from_bus == bu and l.to_bus == bv)
                    or (l.from_bus == bv and l.to_bus == bu)
                ):
                    lines.append(
                        Line(
                            line_id=l.line_id,
                            from_bus=l.from_bus,
                            to_bus=l.to_bus,
                            length_m=l.length_m,
                            cable_type=l.cable_type,
                            r_ohm_per_km=l.r_ohm_per_km,
                            x_ohm_per_km=l.x_ohm_per_km,
                            max_current_a=l.max_current_a,
                        )
                    )
                    break

        grids.append(Grid(nodes=nodes, lines=lines))
        
    grid_system.grids = grids
    return grid_system

def populate_grid_system(
    grid: GridSystem,
    consumption_path: str,
    production_path: str,
    *,
    voltage_lv: float,
    confidence_factor: float,
    cable_catalogue_path: str = Path(__file__).parent / "standard_cables_lv.json",
):
    consumption = load_consumption(consumption_path)
    production = load_production(production_path)
    cables = load_cable_catalogue(cable_catalogue_path)
    for g in grid.grids:
        assign_node_timeseries(g, consumption, production)
        trees = build_directed_trees(g)
        aggregate_downstream_power(g, trees)
        compute_line_currents(g, trees, voltage_lv, confidence_factor)
        assign_cables_to_lines(g, cables)

    return grid


def load_consumption(path) -> dict[int, np.ndarray]:
    df = path if isinstance(path, pd.DataFrame) else pd.read_csv(path, sep=";")
    data = {}

    for col in df.columns:
        try:
            bid = int(col)
        except ValueError:
            continue
        data[bid] = df[col].to_numpy(dtype=float)

    return data


def load_production(path) -> dict[int, np.ndarray]:
    df = path if isinstance(path, pd.DataFrame) else pd.read_csv(path, sep=";")
    data = {}

    for col in df.columns:
        try:
            bid = int(col)
        except ValueError:
            continue
        data[bid] = df[col].to_numpy(dtype=float)

    return data





def load_cable_catalogue(path: str ) -> dict:
    with open(path, "r") as f:
        data = json.load(f)
    return data



def assign_node_timeseries(
    grid: GridSystem,
    consumption: dict[int, np.ndarray],
    production: dict[int, np.ndarray],
):
    for node in grid.nodes:
        if node.node_type == "building":
            cid = node.building_id

            pc = consumption.get(cid)
            pp = production.get(cid)
            
            if pc is None:
                raise ValueError(f"missing consumption for building {cid}")

            if pp is None:
                pp = np.zeros_like(pc)

            node.building_consumption = pc
            node.building_production = pp
            node.consumption = pc
            node.production = pp
            node.net_power = pc - pp

        else:
            node.consumption = None
            node.production = None
            node.net_power = None



def build_directed_trees(
    grid: GridSystem,
) -> dict[int, nx.DiGraph]:
    G = nx.Graph()

    for line in grid.lines.values() if isinstance(grid.lines, dict) else grid.lines:
        u = line.from_bus
        v = line.to_bus
        G.add_edge(u, v)

    trees = {}

    for node in grid.nodes:
        if node.node_type != "supply":
            continue

        root = node.bus_id
        T = nx.DiGraph()

        for parent, child in nx.bfs_edges(G, root):
            T.add_edge(parent, child)

        trees[root] = T

    return trees



def aggregate_downstream_power(
    grid: GridSystem,
    trees: dict[int, nx.DiGraph],
):
    node_by_bus = {n.bus_id: n for n in grid.nodes}
    ref_len = None
    for n in grid.nodes:
        n.downstream_power = 0
        if n.node_type == "building":
            ref_len = n.consumption.shape[0]
    if ref_len is None:
        raise ValueError("no building nodes found")
    zero = np.zeros(ref_len, dtype=float)

    for root, T in trees.items():
        order = list(nx.topological_sort(T))
        order.reverse()

        for u in order:
            node = node_by_bus[u]

            if node.node_type == "building":
                p = node.consumption.astype(float).copy()
            else:
                p = zero.copy()

            for v in T.successors(u):
                child = node_by_bus[v]
                if child.downstream_power is not None:
                    p += child.downstream_power

            node.downstream_power = p



def compute_line_currents(
    grid: GridSystem,
    trees: dict[int, nx.DiGraph],
    voltage_lv: float,
    confidence_factor: float,
):
    sqrt3 = np.sqrt(3.0)
    node_by_bus = {n.bus_id: n for n in grid.nodes}

    for root, T in trees.items():
        for u, v in T.edges:
            parent = node_by_bus[u]
            child = node_by_bus[v]

            p = child.downstream_power.copy()
            p[p < 0.0] = 0.0

            i_ts = (p) / (sqrt3 * voltage_lv)
            i_ts *= confidence_factor
            i_max = float(np.max(i_ts))
            
            for line in grid.lines.values() if isinstance(grid.lines, dict) else grid.lines:
                # line.max_current = 0
                if (
                    (line.from_bus == u and line.to_bus == v)
                    or (line.from_bus == v and line.to_bus == u)
                ):
                    line.current_ts = i_ts
                    line.max_current = i_max
                    break



# def assign_cables_to_lines(
#     grid: GridSystem,
#     cable_catalogue,
#     *,
#     max_parallels: int = 3,
# ):
#     candidates = []
#     for entry in cable_catalogue:
#         i_max = float(entry["max_current_a"])
#         section = float(entry["section_mm2"])
#         r_km = float(entry["r_ohm_per_km"])
#         candidates.append((i_max, section, r_km, entry))

#     lines = grid.lines.values() if isinstance(grid.lines, dict) else grid.lines

#     for line in lines:
#         if not hasattr(line, "max_current"):
#             line.max_current = 0.0

#         i_req = float(line.max_current)

#         best = None
#         best_key = None

#         for i_max, section, r_km, data in candidates:
#             for n in range(1, max_parallels + 1):
#                 if n * i_max < i_req:
#                     continue

#                 r_eff = r_km / n
#                 total_section = n * section

#                 key = (r_eff, total_section)

#                 if best is None or key < best_key:
#                     best = (n, data)
#                     best_key = key

#         if best is None:
#             raise ValueError(
#                 f"no cable configuration can carry {i_req:.2f} A on line {line.line_id}"
#             )

#         n_parallel, data = best

#         line.cable_id = int(data["id"])
#         line.cable_name = data["name"]
#         line.cable_material = data["material"]
#         line.cable_section_mm2 = float(data["section_mm2"])
#         line.cable_max_current = float(data["max_current_a"])
#         line.n_parallel = n_parallel

#         r_per_m = float(data["r_ohm_per_km"]) / 1000.0
#         x_per_m = float(data["x_ohm_per_km"]) / 1000.0

#         line.resistance_per_m = r_per_m / n_parallel
#         line.reactance_per_m = x_per_m / n_parallel

#         line.resistance = line.resistance_per_m * line.length_m
#         line.reactance = line.reactance_per_m * line.length_m

#         line.cable_type = f"{line.cable_name} {line.cable_material}"
#         line.r_ohm_per_km = line.resistance_per_m * 1000.0
#         line.x_ohm_per_km = line.reactance_per_m * 1000.0
#         line.max_current_a = line.cable_max_current * line.n_parallel


def assign_cables_to_lines(
    grid: GridSystem,
    cable_catalogue,
    *,
    max_parallels: int = 3,
):
    candidates = []
    for entry in cable_catalogue:
        i_max = float(entry["max_current_a"])
        section = float(entry["section_mm2"])
        r_km = float(entry["r_ohm_per_km"])
        cost_per_m = float(entry["cost_per_m"])
        candidates.append((i_max, section, r_km, cost_per_m, entry))

    lines = grid.lines.values() if isinstance(grid.lines, dict) else grid.lines

    for line in lines:
        if not hasattr(line, "max_current"):
            line.max_current = 0.0

        i_req = float(line.max_current)

        best = None
        best_cost = None

        for i_max, section, r_km, cost_per_m, data in candidates:
            for n in range(1, max_parallels + 1):
                if n * i_max < i_req:
                    continue

                total_cost = n * cost_per_m * line.length_m

                if best is None or total_cost < best_cost:
                    best = (n, data)
                    best_cost = total_cost

        if best is None:
            raise ValueError(
                f"no cable configuration can carry {i_req:.2f} A on line {line.line_id}"
            )

        n_parallel, data = best

        line.cable_id = int(data["id"])
        line.cable_name = data["name"]
        line.cable_material = data["material"]
        line.cable_section_mm2 = float(data["section_mm2"])
        line.cable_max_current = float(data["max_current_a"])
        line.n_parallel = n_parallel

        r_per_m = float(data["r_ohm_per_km"]) / 1000.0
        x_per_m = float(data["x_ohm_per_km"]) / 1000.0

        line.resistance_per_m = r_per_m / n_parallel
        line.reactance_per_m = x_per_m / n_parallel

        line.resistance = line.resistance_per_m * line.length_m
        line.reactance = line.reactance_per_m * line.length_m

        line.cable_type = f"{line.cable_name} {line.cable_material}"
        line.r_ohm_per_km = line.resistance_per_m * 1000.0
        line.x_ohm_per_km = line.reactance_per_m * 1000.0
        line.max_current_a = line.cable_max_current * line.n_parallel
