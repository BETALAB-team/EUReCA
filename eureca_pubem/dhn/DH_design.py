import geopandas as gpd
import numpy as np
import networkx as nx
import pandas as pd 
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import split, unary_union, nearest_points

from eureca_pubem.dhn.DH_objects import Node, Line

def load_dhn_demand(csv_path):
    import pandas as pd
    df = csv_path if isinstance(csv_path, pd.DataFrame) else pd.read_csv(csv_path, sep=";")
    return df
def get_buildings_to_attach(scenario, peak_demands):
    demanded = {bid for bid, v in peak_demands.items() if v > 0.0}

    served = {
        n.node_id
        for dhn in scenario.District_Heating_Systems.values()
        for n in dhn.nodes
        if n.node_type == "consumer" and n.active
    }

    return demanded - served
def find_closest_dhn_element(dhn, point):
    import numpy as np

    def dist(p, q):
        return np.hypot(p[0] - q[0], p[1] - q[1])

    def project(p, a, b):
        ap = np.array(p) - np.array(a)
        ab = np.array(b) - np.array(a)
        t = np.dot(ap, ab) / np.dot(ab, ab)
        t = max(0.0, min(1.0, t))
        proj = np.array(a) + t * ab
        return proj, np.linalg.norm(ap - t * ab)

    nodes = {n.node_id: n for n in dhn.nodes}
    active_nodes = [n for n in dhn.nodes if n.active]
    active_lines = [l for l in dhn.lines if l.active]

    best = None
    best_dist = np.inf

    for n in active_nodes:
        d = dist(point, (n.x, n.y))
        if d < best_dist:
            best_dist = d
            best = ("node", n.node_id)

    for l in active_lines:
        a = nodes[l.start_node]
        b = nodes[l.end_node]
        proj, d = project(point, (a.x, a.y), (b.x, b.y))
        if d < best_dist:
            best_dist = d
            best = ("line", l, proj)

    return best, best_dist

def split_line(dhn, line, x, y):
    nodes = {n.node_id: n for n in dhn.nodes}

    next_node_id = max(nodes) + 1
    next_line_id = max(l.line_id for l in dhn.lines) + 1

    origin = getattr(line, "origin_line_id", line.line_id)

    # new connection node
    cn = type(dhn.nodes[0])(
        node_id=next_node_id,
        node_type="connection",
    )
    cn.x = x
    cn.y = y
    cn.active = True
    cn.is_new = True
    dhn.nodes.append(cn)

    # deactivate old line
    line.active = False

    # A–C
    l1 = type(line)(
        line_id=next_line_id,
        start_node=line.start_node,
        end_node=next_node_id,
        length=((nodes[line.start_node].x - x)**2 + (nodes[line.start_node].y - y)**2) ** 0.5,
    )
    l1.active = True
    l1.is_new = True
    l1.origin_line_id = origin
    dhn.lines.append(l1)
    next_line_id += 1

    # C–B
    l2 = type(line)(
        line_id=next_line_id,
        start_node=next_node_id,
        end_node=line.end_node,
        length=((nodes[line.end_node].x - x)**2 + (nodes[line.end_node].y - y)**2) ** 0.5,
    )
    l2.active = True
    l2.is_new = True
    l2.origin_line_id = origin
    dhn.lines.append(l2)

    return next_node_id

def attach_single_building(dhn, building, demand_array, closest, distance):
    nodes = {n.node_id: n for n in dhn.nodes}
    next_line_id = max(l.line_id for l in dhn.lines) + 1

    if closest[0] == "node":
        target_id = closest[1]
    else:
        line, (x, y) = closest[1], closest[2]
        target_id = split_line(dhn, line, x, y)

    # new consumer node
    bn = type(dhn.nodes[0])(
        node_id=building.building_id,
        node_type="consumer",
    )
    bn.x = building.x
    bn.y = building.y
    bn.demand = demand_array
    bn.active = True
    bn.is_new = True
    dhn.nodes.append(bn)

    # service line
    sl = type(dhn.lines[0])(
        line_id=next_line_id,
        start_node=target_id,
        end_node=bn.node_id,
        length=distance,
    )
    sl.active = True
    sl.is_new = True
    sl.origin_line_id = None  # service pipe
    dhn.lines.append(sl)

def attach_missing_dhn_consumers(
    scenario,
    peak_demands: dict[int, float],
    csv_path,
):
    df = load_dhn_demand(csv_path)
    to_attach = get_buildings_to_attach(scenario, peak_demands)

    for bid in to_attach:
        b = scenario.buildings.get(bid)
        if b is None:
            continue

        best = None
        best_dist = float("inf")
        best_dhn = None

        # IMPORTANT: search on CURRENT tree
        for dhn in scenario.District_Heating_Systems.values():
            res, d = find_closest_dhn_element(dhn, (b.x, b.y))
            if res is not None and d < best_dist:
                best = res
                best_dist = d
                best_dhn = dhn

        if best_dhn is None:
            continue

        attach_single_building(
            best_dhn,
            b,
            df[str(bid)].to_numpy(),
            best,
            best_dist,
        )

        # CRITICAL: rebuild tree after EACH attachment
        rebuild_spanning_tree_index(best_dhn)

def rebuild_spanning_tree_index(system):
    from collections import deque

    system.spanning_tree_index = []

    nodes = system.nodes
    node_ids = [n.node_id for n in nodes]
    id_to_pos = {nid: i for i, nid in enumerate(node_ids)}

    supply = [n.node_id for n in nodes if n.node_type == "supply"]
    if len(supply) != 1:
        raise RuntimeError("Expected exactly one supply node")

    root = id_to_pos[supply[0]]

    adj = [[] for _ in nodes]
    for eidx, l in enumerate(system.lines):
        if not l.active:
            continue
        u = id_to_pos[l.start_node]
        v = id_to_pos[l.end_node]
        adj[u].append((v, eidx))
        adj[v].append((u, eidx))

    visited = [False] * len(nodes)
    visited[root] = True

    q = deque([root])
    tree_edges = []

    while q:
        u = q.popleft()
        for v, eidx in adj[u]:
            if visited[v]:
                continue
            visited[v] = True
            tree_edges.append(eidx)
            q.append(v)

    reachable = sum(visited)

    if len(tree_edges) != reachable - 1:
        raise RuntimeError(
            "Topology is not a tree after attachment "
            f"(edges={len(tree_edges)}, nodes={reachable})"
        )

    system.spanning_tree_index = tree_edges


def size_all_dhn_systems(scenario):
    from dhn.DH_design import lines_sizing, rebuild_spanning_tree_index

    for system in scenario.District_Heating_Systems.values():
        rebuild_spanning_tree_index(system)
        lines_sizing(system)


def read_peak_demands(csv_path):
    import pandas as pd
    df = csv_path if isinstance(csv_path, pd.DataFrame) else pd.read_csv(csv_path, sep=";")

    peaks = {}
    for col in df.columns:
        try:
            bid = int(col)
        except Exception:
            continue
        peaks[bid] = float(df[col].max())

    return peaks

def apply_dhn_activation_from_demands(
    dhn,
    peak_demands: dict[int, float],
):
    from collections import deque

    for n in dhn.nodes:
        n.active = True
        n.downstream_demand = 0.0

    for l in dhn.lines:
        l.active = True

    nodes = {n.node_id: n for n in dhn.nodes}

    adj = {nid: [] for nid in nodes}
    for l in dhn.lines:
        if l.start_node in nodes and l.end_node in nodes:
            adj[l.start_node].append(l.end_node)
            adj[l.end_node].append(l.start_node)

    supplies = [n.node_id for n in dhn.nodes if n.node_type == "supply"]
    if not supplies:
        return

    parent = {}
    children = {nid: [] for nid in nodes}

    q = deque(supplies)
    visited = set(supplies)

    while q:
        u = q.popleft()
        for v in adj[u]:
            if v in visited:
                continue
            visited.add(v)
            parent[v] = u
            children[u].append(v)
            q.append(v)

    order = []

    def dfs(u):
        for v in children[u]:
            dfs(v)
        order.append(u)

    for s in supplies:
        dfs(s)

    for nid in order:
        n = nodes[nid]
        if n.node_type == "consumer":
            n.downstream_demand = peak_demands.get(nid, 0.0)
        else:
            n.downstream_demand = sum(
                nodes[c].downstream_demand for c in children[nid]
            )

    for l in dhn.lines:
        if l.end_node in nodes:
            if nodes[l.end_node].downstream_demand == 0.0:
                l.active = False
                nodes[l.end_node].active = False

    for s in supplies:
        nodes[s].active = True
        





def parse_buildings_nodes(buildings_path, streets_path, 
                          street_flag = "path",
                          building_flag = "path"):
    def load_gdf(obj):
        if isinstance(obj, (str, bytes)):
            return gpd.read_file(obj)
        elif isinstance(obj, gpd.GeoDataFrame):
            return obj
        else:
            raise TypeError("Input must be a GeoDataFrame or a path to a file.") 
    streets = load_gdf (streets_path)

    def load_gdf(obj):
        if isinstance(obj, (str, bytes)):
            return gpd.read_file(obj)
        elif isinstance(obj, gpd.GeoDataFrame):
            return obj
        else:
            raise TypeError("Input must be a GeoDataFrame or a path to a file.")    

    buildings = load_gdf(buildings_path)


    if streets.crs is None:
        streets = streets.set_crs(buildings.crs)
    if buildings.crs is None:
        buildings = buildings.set_crs(streets.crs)

    ref = streets.to_crs(4326) if streets.crs.to_epsg() != 4326 else streets
    c = ref.geometry.unary_union.centroid
    zone = int((c.x + 180) // 6) + 1
    crs = f"EPSG:{32600 + zone}"

    streets = streets.to_crs(crs)
    buildings = buildings.to_crs(crs)

    nodes = {}

    building_ids = set(int(row["id"]) for _, row in buildings.iterrows())
    max_building_id = max(building_ids)
    next_connection_id = max_building_id + 1

    supply_nid = None

    for _, row in buildings.iterrows():
        bid = int(row["id"])
        t = int(row["type"])
        if t not in (-1, 1):
            continue
        ntype = "supply" if t == 1 else "consumer"
        n = Node(node_id=bid, node_type=ntype)
        n.active = True
        n.x = row.geometry.x
        n.y = row.geometry.y
        if ntype == "supply" and "capacity" in row:
            n.capacity = float(row["capacity"])
        nodes[bid] = n
        if t == 1:
            supply_nid = bid
        
    return nodes, next_connection_id


def build_dhn_from_geometry(
    buildings_path,
    streets_path,
    *,
    demand_csv_path,
    crs_metric=None,
    grid_resolution=None,
    grid_penalty=None,
):

    buildings_gdf, streets_gdf = load_and_project_inputs(
        buildings_path,
        streets_path,
        crs_metric,
    )

    terminals = parse_buildings(buildings_gdf)

    peak_demands = load_peak_demands_from_csv(demand_csv_path)
    attach_peak_demands(terminals, peak_demands)

    streets = sanitize_street_geometry(streets_gdf)

    street_graph = build_street_graph(streets)
    aux_graph = build_auxiliary_grid(streets, grid_resolution, grid_penalty)
    routing = merge_graphs(street_graph, aux_graph)
    routing = attach_terminal_nodes(routing, terminals)

    assignment = assign_buildings_to_supply(
        routing,
        terminals,
    )

    trees = extract_supply_trees(
        routing,
        terminals,
        assignment,
    )
    nodes, lines = flatten_trees_to_pruned_output(
        trees,
        terminals,
    )

    return nodes, lines

def load_and_project_inputs(
    buildings_path,
    streets_path,
    crs_metric=None,
):
    buildings_gdf = (
        gpd.read_file(buildings_path)
        if not isinstance(buildings_path, gpd.GeoDataFrame)
        else buildings_path.copy()
    )
    streets_gdf = (
        gpd.read_file(streets_path)
        if not isinstance(streets_path, gpd.GeoDataFrame)
        else streets_path.copy()
    )

    if buildings_gdf.crs is None or streets_gdf.crs is None:
        raise ValueError("missing CRS")

    if crs_metric is None:
        ref = buildings_gdf.to_crs(4326)
        c = ref.geometry.unary_union.centroid
        zone = int((c.x + 180) // 6) + 1
        crs_metric = f"EPSG:{32600 + zone}"

    buildings_gdf = buildings_gdf.to_crs(crs_metric)
    streets_gdf = streets_gdf.to_crs(crs_metric)

    return buildings_gdf, streets_gdf

def parse_buildings(buildings_gdf):
    terminals = {}
    for _, row in buildings_gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        nid = int(row["id"])
        t = int(row["type"])
        if t not in (-1, 1):
            continue

        node_type = "supply" if t == 1 else "consumer"

        n = Node(node_id=nid, node_type=node_type)
        n.active = True
        n.x = float(geom.x)
        n.y = float(geom.y)

        if node_type == "supply" and "capacity" in row:
            try:
                n.capacity = float(row["capacity"])
            except Exception:
                pass

        terminals[nid] = n

    return terminals

def load_peak_demands_from_csv(path):
    import pandas as pd

    df = path if isinstance(path, pd.DataFrame) else pd.read_csv(path)

    peaks = {}
    for col in df.columns:
        try:
            bid = int(col)
        except Exception:
            continue
        peaks[bid] = float(df[col].max())
    return peaks


def attach_peak_demands(
    terminals: dict[int, Node],
    peak_demands: dict[int, float],
):
    for nid, n in terminals.items():
        if n.node_type != "consumer":
            continue

        val = peak_demands.get(nid, 0.0)
        try:
            val = float(np.asarray(val).reshape(-1)[0])
        except Exception:
            val = 0.0

        if not np.isfinite(val) or val < 0.0:
            val = 0.0
            
        n.design_demand = val



def sanitize_street_geometry(
    streets_gdf,
    *,
    tol: float = 1e-6,
):
    if not isinstance(streets_gdf, gpd.GeoDataFrame):
        streets_gdf = gpd.GeoDataFrame(
            streets_gdf,
            geometry="geometry",
            crs=getattr(streets_gdf, "crs", None),
        )

    crs = streets_gdf.crs
    gdf = streets_gdf.copy()

    gdf["geometry"] = gdf.geometry.apply(
        lambda g: list(g.geoms) if isinstance(g, MultiLineString) else [g]
    )
    gdf = gdf.explode("geometry", ignore_index=True)
    gdf = gdf[gdf.geometry.apply(lambda g: g is not None and not g.is_empty)]
    gdf = gdf[gdf.geometry.apply(lambda g: g.length > tol)].reset_index(drop=True)

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
                pass

        if isinstance(geom, MultiLineString):
            split_lines.extend(list(geom.geoms))
        else:
            split_lines.append(geom)

    out = gpd.GeoDataFrame(
        geometry=[g for g in split_lines if g.length > tol],
        crs=crs,
    ).reset_index(drop=True)

    return out

def build_street_graph(streets_gdf):
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

                w = float(np.hypot(v[0] - u[0], v[1] - u[1]))

                if G.has_edge(u, v):
                    if w < G[u][v]["weight"]:
                        G[u][v]["weight"] = w
                else:
                    G.add_edge(u, v, weight=w, kind="street")

    return G


def build_auxiliary_grid(
    streets_gdf,
    grid_resolution=None,
    grid_penalty=None,
):
    xmin, ymin, xmax, ymax = streets_gdf.total_bounds

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
                vx = x + dx
                vy = y + dy
                if vx < xmin or vx > xmax or vy < ymin or vy > ymax:
                    continue
                v = (float(vx), float(vy))
                if v not in G:
                    continue
                w = grid_penalty * float(np.hypot(dx, dy))
                G.add_edge(u, v, weight=w, kind="grid")

    return G


def merge_graphs(
    street_graph: nx.Graph,
    aux_graph: nx.Graph,
):
    G = nx.Graph()

    for u, data in street_graph.nodes(data=True):
        G.add_node(u, **data)

    for u, v, data in street_graph.edges(data=True):
        G.add_edge(u, v, **data)

    for u, data in aux_graph.nodes(data=True):
        if u not in G:
            G.add_node(u, **data)

    for u, v, data in aux_graph.edges(data=True):
        if G.has_edge(u, v):
            if data.get("weight", np.inf) < G[u][v].get("weight", np.inf):
                G[u][v].update(data)
        else:
            G.add_edge(u, v, **data)

    street_nodes = np.array(list(street_graph.nodes), dtype=float)
    grid_nodes = np.array(list(aux_graph.nodes), dtype=float)

    for sn in street_graph.nodes:
        d = grid_nodes - np.array(sn)
        dist = np.hypot(d[:, 0], d[:, 1])
        i = int(np.argmin(dist))
        gn = tuple(grid_nodes[i])

        if not G.has_edge(sn, gn):
            G.add_edge(sn, gn, weight=float(dist[i]), kind="street_grid_link")

    return G

def attach_terminal_nodes(
    graph: nx.Graph,
    terminals: dict[int, Node],
):
    G = graph.copy()

    graph_nodes = np.array(list(G.nodes), dtype=float)

    for n in terminals.values():
        p = np.array([n.x, n.y], dtype=float)

        d = graph_nodes - p
        dist = np.hypot(d[:, 0], d[:, 1])
        i = int(np.argmin(dist))

        nearest = tuple(graph_nodes[i])
        u = (float(p[0]), float(p[1]))

        if u not in G:
            G.add_node(u)

        G.add_edge(u, nearest, weight=float(dist[i]), kind="terminal")

    return G

def assign_buildings_to_supply(
    graph: nx.Graph,
    terminals: dict[int, Node],
):
    supplies = {
        nid: n for nid, n in terminals.items()
        if n.node_type == "supply"
    }
    consumers = {
        nid: n for nid, n in terminals.items()
        if n.node_type == "consumer"
    }

    remaining_capacity = {}
    for sid, s in supplies.items():
        cap = getattr(s, "capacity", None)
        remaining_capacity[sid] = float(cap) if cap is not None else np.inf

    assignment = {}

    for cid, c in consumers.items():
        c_coord = (c.x, c.y)

        best_sid = None
        best_cost = np.inf

        demand = getattr(c, "design_demand", 0.0)
        if demand is None:
            demand = 0.0

        for sid, s in supplies.items():
            if remaining_capacity[sid] < demand:
                continue

            s_coord = (s.x, s.y)

            try:
                d = nx.shortest_path_length(
                    graph,
                    source=s_coord,
                    target=c_coord,
                    weight="weight",
                )
            except nx.NetworkXNoPath:
                continue

            if d < best_cost:
                best_cost = d
                best_sid = sid

        if best_sid is None:
            continue

        assignment[cid] = best_sid
        remaining_capacity[best_sid] -= demand
        
    assigned_ids = set(assignment.keys())

    for cid, c in consumers.items():
        if cid not in assigned_ids:
            c.active = False

    return assignment

def extract_supply_trees(
    graph: nx.Graph,
    terminals: dict[int, Node],
    assignment: dict[int, int],
):
    supplies = {
        nid: n for nid, n in terminals.items()
        if n.node_type == "supply"
    }

    trees = {}

    for sid, s in supplies.items():
        s_coord = (s.x, s.y)
        Gs = nx.Graph()

        consumers = [
            cid for cid, asid in assignment.items()
            if asid == sid
        ]

        for cid in consumers:
            c = terminals[cid]
            c_coord = (c.x, c.y)

            try:
                path = nx.shortest_path(
                    graph,
                    source=s_coord,
                    target=c_coord,
                    weight="weight",
                )
            except nx.NetworkXNoPath:
                continue

            for u, v in zip(path[:-1], path[1:]):
                w = graph[u][v]["weight"]
                kind = graph[u][v].get("kind", None)
                if not Gs.has_edge(u, v):
                    Gs.add_edge(u, v, weight=w, kind=kind)

        if Gs.number_of_nodes() == 0:
            trees[sid] = nx.Graph()
            continue

        T = nx.minimum_spanning_tree(Gs, weight="weight")
        trees[sid] = T

    return trees

def flatten_trees_to_pruned_output(
    trees: dict[int, nx.Graph],
    terminals: dict[int, Node],
):
    used_coords = set()
    edge_set = set()

    for T in trees.values():
        for u, v in T.edges:
            edge_set.add((u, v))
            edge_set.add((v, u))
            used_coords.add(u)
            used_coords.add(v)

    coord_to_node = {}
    pruned_nodes = []

    for n in terminals.values():
        coord = (float(n.x), float(n.y))
        coord_to_node[coord] = n
        pruned_nodes.append(n)
        used_coords.add(coord)

    next_id = max(n.node_id for n in pruned_nodes) + 1 if pruned_nodes else 1

    for coord in used_coords:
        if coord in coord_to_node:
            continue
        n = Node(node_id=next_id, node_type="connection")
        n.active = True
        n.x = coord[0]
        n.y = coord[1]
        coord_to_node[coord] = n
        pruned_nodes.append(n)
        next_id += 1

    pruned_lines = []
    lid = 1

    for u, v in edge_set:
        nu = coord_to_node[u].node_id
        nv = coord_to_node[v].node_id
        if nu == nv:
            continue
        length = float(np.hypot(v[0] - u[0], v[1] - u[1]))
        l = Line(
            line_id=lid,
            start_node=nu,
            end_node=nv,
            length=length,
        )
        l.active = True
        pruned_lines.append(l)

        lid += 1

    return pruned_nodes, pruned_lines
