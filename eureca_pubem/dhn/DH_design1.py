import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt
import math

from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import unary_union, split, nearest_points

from dhn.DH_objects import Node, Line


def parse_buildings_nodes(buildings_path, streets_path, 
                          street_flag = "path",
                          building_flag = "path"):
    if street_flag == "path":
        streets = gpd.read_file(streets_path)
    else:
        streets = streets_path
        
    if building_flag == "path":    
        buildings = gpd.read_file(buildings_path)
    else:
        buildings = buildings_path

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
        n.x = row.geometry.x
        n.y = row.geometry.y
        if ntype == "supply" and "capacity" in row:
            n.capacity = float(row["capacity"])
        nodes[bid] = n
        if t == 1:
            supply_nid = bid
        
    return nodes, next_connection_id

def build_prune_plot_dh(
    streets_path,
    buildings_path,
    alpha=2e-5,
    street_flag = "path",
    building_flag = "path"
):
    if street_flag == "path":
        streets = gpd.read_file(streets_path)
    else:
        streets = streets_path
        
    if building_flag == "path":    
        buildings = gpd.read_file(buildings_path)
    else:
        buildings = buildings_path

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
    lines = []

    # building_ids = set(int(row["id"]) for _, row in buildings.iterrows())
    # max_building_id = max(building_ids)
    # next_connection_id = max_building_id + 1

    # supply_nid = None

    # for _, row in buildings.iterrows():
    #     bid = int(row["id"])
    #     t = int(row["type"])
    #     if t not in (-1, 1):
    #         continue
    #     ntype = "supply" if t == 1 else "consumer"
    #     n = Node(node_id=bid, node_type=ntype)
    #     n.x = row.geometry.x
    #     n.y = row.geometry.y
    #     if ntype == "supply" and "capacity" in row:
    #         n.capacity = float(row["capacity"])
    #     nodes[bid] = n
    #     if t == 1:
    #         supply_nid = bid
    nodes, next_connection_id = parse_buildings_nodes(buildings_path, streets_path, building_flag = "direct_dataframe")
    

    for ident, n in nodes.items():
        bid = ident
        t = n.node_type
        n = Node(node_id=bid, node_type=t)
        if t == "supply":
            supply_nid = bid
   
            
    def new_connection_node(x, y):
        nonlocal next_connection_id
        nid = next_connection_id
        n = Node(node_id=nid, node_type="connection")
        n.x = x
        n.y = y
        nodes[nid] = n
        next_connection_id += 1
        return nid

    def iter_lines(g):
        if isinstance(g, LineString):
            yield g
        elif isinstance(g, MultiLineString):
            for x in g.geoms:
                yield x

    all_lines = []
    for g in streets.geometry:
        for l in iter_lines(g):
            all_lines.append(l)

    merged = unary_union(all_lines)

    split_lines = []
    for l in all_lines:
        try:
            parts = split(l, merged)
            split_lines.extend(list(parts.geoms))
        except Exception:
            split_lines.append(l)

    street_node_map = {}
    street_segments = []
    lid = 1

    def street_node(p):
        key = (round(p.x, 3), round(p.y, 3))
        if key not in street_node_map:
            street_node_map[key] = new_connection_node(p.x, p.y)
        return street_node_map[key]

    for idx, l in enumerate(split_lines):
        a = Point(l.coords[0])
        b = Point(l.coords[-1])
        na = street_node(a)
        nb = street_node(b)

        w = streets.iloc[min(idx, len(streets) - 1)].get("weight", 1.0)
        street_segments.append((l, w))

        lines.append(
            Line(
                line_id=lid,
                start_node=na,
                end_node=nb,
                length=a.distance(b),
            )
        )
        lid += 1

    for _, row in buildings.iterrows():
        t = int(row["type"])
        if t not in (-1, 1):
            continue

        bid = int(row["id"])
        geom = row.geometry

        seg, seg_w = min(street_segments, key=lambda s: geom.distance(s[0]))
        p = nearest_points(geom, seg)[1]

        cn = new_connection_node(p.x, p.y)

        lines.append(
            Line(
                line_id=lid,
                start_node=bid,
                end_node=cn,
                length=geom.distance(p),
            )
        )
        lid += 1

        a = Point(seg.coords[0])
        b = Point(seg.coords[-1])
        na = street_node(a)
        nb = street_node(b)

        lines.append(Line(lid, na, cn, a.distance(p)))
        lid += 1
        lines.append(Line(lid, cn, nb, b.distance(p)))
        lid += 1

    G = nx.Graph()
    supply_point = Point(nodes[supply_nid].x, nodes[supply_nid].y)

    for l in lines:
        n1 = nodes[l.start_node]
        n2 = nodes[l.end_node]
        mid = Point((n1.x + n2.x) / 2, (n1.y + n2.y) / 2)
        d = mid.distance(supply_point)
        w = l.length * math.exp(alpha * d)
        G.add_edge(l.start_node, l.end_node, weight=w)

    useful = set()

    for n in nodes.values():
        if n.node_type == "consumer":
            try:
                path = nx.shortest_path(G, supply_nid, n.node_id, weight="weight")
            except nx.NetworkXNoPath:
                continue
            for u, v in zip(path[:-1], path[1:]):
                useful.add((u, v))
                useful.add((v, u))

    pruned_lines = []
    used_nodes = set()

    for l in lines:
        if (l.start_node, l.end_node) in useful:
            pruned_lines.append(l)
            used_nodes.add(l.start_node)
            used_nodes.add(l.end_node)

    pruned_nodes = [
        n for n in nodes.values()
        if n.node_type in ("consumer", "supply") or n.node_id in used_nodes
    ]

    return pruned_nodes, pruned_lines


