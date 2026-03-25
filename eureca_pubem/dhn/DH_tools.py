import numpy as np
import pandas as pd
from collections import deque
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec
from eureca_pubem.dhn.pipe_catalog_loader import (
    DN_CATALOG,
    DN_INNER_D,
    DN_OUTER_D,
    PIPE_MATERIALS,
    INSULATION_CLASSES,
)
try:
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False
    


import math

def size_dhn_with_baseline(
    baseline,
    scenario,
    **sizing_kwargs,
):
    """
    Sizes DHN lines in `scenario` using `baseline` as reference.

    Rules:
    - new line            -> keep sized value
    - existing line:
        - if new size > baseline size -> upgrade
        - else -> keep baseline size

    Output stored in:
        scenario.dhn_pipe_changes = {
            "new_pipes": [...],
            "pipes_changed": [(old_line, new_line), ...],
        }
    """

    # run sizing on BOTH
    # for sys_id, sys_base in baseline.District_Heating_Systems.items():
    #     lines_sizing(sys_base, **sizing_kwargs)

    for sys_id, sys_new in scenario.District_Heating_Systems.items():
        rebuild_spanning_tree_index(sys_new)
        lines_sizing(sys_new, **sizing_kwargs)

    new_pipes = []
    pipes_changed = []

    def line_key(l):
        return tuple(sorted((l.start_node, l.end_node)))

    for sys_id, sys_new in scenario.District_Heating_Systems.items():
        sys_base = baseline.District_Heating_Systems[sys_id]

        base_lines = {
            line_key(l): l
            for l in sys_base.lines
        }

        for l_new in sys_new.lines:
            k = line_key(l_new)
            if getattr(l_new, "is_new", False) or k not in base_lines:
                new_pipes.append(l_new)
                continue

            l_old = base_lines[k]

            if l_new.inner_diameter_required > l_old.inner_diameter:
                pipes_changed.append((l_old, l_new))
            else:
                l_new = l_old

    scenario.dhn_pipe_changes = {
        "new_pipes": new_pipes,
        "pipes_changed": pipes_changed,
    }


def assign_consumers_to_supplies(
    nodes: dict,
    *,
    concurrency_fn,
    max_distance_fn=None,
    huge_distance: float = 1e9,
    fail_on_unassigned: bool = False,
):
    """
    Assign consumer nodes to supply nodes based on proximity and capacity.

    Mutates:
        node.assigned_supply_id

    Returns:
        dict[supply_id, list[consumer_ids]]
    """

    supplies = []
    consumers = []

    for n in nodes.values():
        if n.node_type == "supply":
            if not hasattr(n, "capacity") or not math.isfinite(n.capacity):
                raise ValueError(f"Supply node {n.node_id} has no valid capacity")
            supplies.append(n)
        elif n.node_type == "consumer":
            if not hasattr(n, "peak_demand") or n.peak_demand < 0 or not math.isfinite(n.peak_demand):
                raise ValueError(f"Consumer node {n.node_id} has no valid peak_demand")

            consumers.append(n)

    if not supplies:
        raise ValueError("No supply nodes found")

    if not consumers:

        out = {}
        for s in supplies:
            s.assigned_supply_id = s.node_id
            out[s.node_id] = []
        return out



    supply_state = {}

    for s in supplies:
        s.assigned_supply_id = s.node_id
        max_dist = (
            max_distance_fn(s)
            if max_distance_fn is not None
            else huge_distance
        )
        supply_state[s.node_id] = {
            "node": s,
            "sum_peak": 0.0,
            "n": 0,
            "max_distance": float(max_dist),
            "consumers": [],
        }


    consumers_sorted = sorted(
        consumers,
        key=lambda c: c.peak_demand,
        reverse=True,
    )



    def dist(a, b):
        dx = float(a.x) - float(b.x)
        dy = float(a.y) - float(b.y)
        return math.hypot(dx, dy)


    unassigned = []

    for c in consumers_sorted:
        c.assigned_supply_id = None


        ordered_supplies = sorted(
            supply_state.values(),
            key=lambda s: dist(c, s["node"]),
        )

        assigned = False

        for s in ordered_supplies:
            d = dist(c, s["node"])
            if d > s["max_distance"]:
                continue

            new_sum = s["sum_peak"] + c.peak_demand
            new_n = s["n"] + 1
            corrected = new_sum * concurrency_fn(new_n)

            if corrected <= s["node"].capacity:
                # assign
                s["sum_peak"] = new_sum
                s["n"] = new_n
                s["consumers"].append(c.node_id)
                c.assigned_supply_id = s["node"].node_id
                assigned = True
                break

        if not assigned:
            unassigned.append(c.node_id)


    if unassigned and fail_on_unassigned:
        raise ValueError(
            f"Unassigned consumers: {unassigned}. "
            "Check capacities, distances, or concurrency."
        )


    result = {}
    for sid, s in supply_state.items():
        result[sid] = list(s["consumers"])

    return result



def incidence_matrix(system) -> np.ndarray:
    node_ids = sorted(n.node_id for n in system.nodes if n.active)
    n_nodes = len(node_ids)
    n_lines = len(system.lines)

    node_index = {nid: i for i, nid in enumerate(node_ids)}

    A = np.zeros((n_nodes, n_lines), dtype=int)

    for j, line in enumerate(system.lines):
        A[node_index[line.start_node], j] = 1
        A[node_index[line.end_node], j] = -1

    connections = np.array(
        [[line.start_node, line.end_node] for line in system.lines if line.active],
        dtype=int
    )
    
    num_nodes = A.shape[0]
    
    tree_idx, link_idx = _spanning_tree_indices(connections)
    A_T, A_L = _split_incidence_matrix(A, tree_idx, link_idx)

    system.incidence_matrix = A
    system.spanning_tree_incidence_matrix = A_T
    system.branches_incidence_matrix = A_L
    system.spanning_tree_index = tree_idx
    system.branches_index = link_idx
    



def _spanning_tree_indices(connections):
    node_ids = set(connections.flatten())
    parent = {int(n): int(n) for n in node_ids}

    def find(x):
        x = int(x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        parent[rb] = ra
        return True

    tree_idx = []
    link_idx = []

    for i, (n1, n2) in enumerate(connections):
        if union(n1, n2):
            tree_idx.append(i)
        else:
            link_idx.append(i)

    return tree_idx, link_idx


def _split_incidence_matrix(A, tree_idx, link_idx):
    A_T = A[:, tree_idx]
    A_L = A[:, link_idx]
    return A_T, A_L


def _rng(seed=None):
    return np.random.default_rng(seed)


def _choose_supply_datum_node(dh, seed=None):
    rng = _rng(seed)

    supply_nodes = [n for n in dh.nodes if getattr(n, "node_type", None) == "supply"]
    caps = []
    has_any_cap = False
    for n in supply_nodes:
        cap = getattr(n, "max_capacity", None)
        if cap is None:
            cap = getattr(n, "capacity", None)
        if cap is None:
            cap = getattr(n, "max_capacity_kw", None)
        if cap is None:
            cap = getattr(n, "max_capacity_mw", None)
        if cap is None:
            caps.append(np.nan)
        else:
            has_any_cap = True
            caps.append(float(cap))

    if not has_any_cap:
        return int(rng.choice([n.node_id for n in supply_nodes]))

    max_cap = np.nanmax(caps)
    top = [supply_nodes[i].node_id for i, c in enumerate(caps) if np.isfinite(c) and c == max_cap]
    if not top:
        return int(rng.choice([n.node_id for n in supply_nodes]))
    return int(rng.choice(top))


def cutset_matrix_with_supply_datum(system, seed=None):
    
    A_T = system.spanning_tree_incidence_matrix
    A_L = system.branches_incidence_matrix

    node_ids = system.node_ids            # MUST be the same used in incidence
    node_index = {nid: i for i, nid in enumerate(node_ids)}

    datum_node_id = _choose_supply_datum_node(system, seed=seed)
    datum_row = node_index[datum_node_id]

    n_nodes = len([n for n in system.nodes if n.active])

    keep = np.array([i for i in range(n_nodes) if i != datum_row], dtype=int)

    A_Tr = A_T[keep, :]
    A_Lr = A_L[keep, :]

    if A_Tr.shape[0] != A_Tr.shape[1]:
        import networkx as nx
        
        G = nx.Graph()
        
        for line in system.lines:
            if line.active:
                G.add_edge(line.start_node, line.end_node)
        

        raise ValueError(
            f"Invalid spanning tree: A_Tr is {A_Tr.shape}, "
            f"expected ({n_nodes-1}, {n_nodes-1})"
        )

    D_L = np.linalg.solve(A_Tr, A_Lr)
    D = np.hstack([np.eye(n_nodes - 1), D_L])

    system.cutset_matrix = D
    system.cutset_links_block = D_L



def loop_matrix(system):
    D_L = system.cutset_links_block  

    B_T = -D_L.T                     
    n_L = B_T.shape[0]

    B = np.hstack([B_T, np.eye(n_L)])

    system.loop_matrix = B
    system.loop_tree_block = B_T


def _value_at_t(x, t):
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        return float(a)
    a = a.reshape(-1)
    return float(a[t])


def build_tree_pressure_cache(system, seed=None, T=None):
    nodes = system.nodes
    lines = system.lines

    node_ids = np.array([int(n.node_id) for n in nodes], dtype=int)
    id_to_pos = {nid: i for i, nid in enumerate(node_ids)}
    n_nodes = len(nodes)

    supply_ids = [int(n.node_id) for n in nodes if getattr(n, "node_type", None) == "supply"]


    rng = np.random.default_rng(seed)
    datum_id = int(getattr(system, "datum_node", rng.choice(supply_ids)))

    datum_pos = id_to_pos[datum_id]

    tree_idx = np.asarray(system.spanning_tree_index, dtype=int)

    adj = [[] for _ in range(n_nodes)]
    for eidx in tree_idx:
        l = lines[int(eidx)]
        u_id = int(l.start_node)
        v_id = int(l.end_node)
        u = id_to_pos[u_id]
        v = id_to_pos[v_id]
        adj[u].append((v, int(eidx)))
        adj[v].append((u, int(eidx)))

    parent = np.full(n_nodes, -1, dtype=int)
    parent_edge = np.full(n_nodes, -1, dtype=int)
    children = [[] for _ in range(n_nodes)]
    order = []

    parent[datum_pos] = datum_pos
    q = deque([datum_pos])
    while q:
        u = q.popleft()
        order.append(u)
        for v, eidx in adj[u]:
            if parent[v] != -1:
                continue
            parent[v] = u
            parent_edge[v] = eidx
            children[u].append(v)
            q.append(v)

    sign = np.zeros(n_nodes, dtype=float)
    for v in range(n_nodes):
        if v == datum_pos:
            continue
        eidx = int(parent_edge[v])
        l = lines[eidx]
        p = int(parent[v])
        v_id = int(node_ids[v])
        p_id = int(node_ids[p])

        if int(l.start_node) == p_id and int(l.end_node) == v_id:
            sign[v] = +1.0
        elif int(l.start_node) == v_id and int(l.end_node) == p_id:
            sign[v] = -1.0

    L = np.array([float(l.length) for l in lines], dtype=float)
    D = np.array([float(l.inner_diameter) for l in lines], dtype=float)/1000
    eps = np.array([float(l.roughness) for l in lines], dtype=float)

    if T is None:
        T = None
        for n in nodes:
            d = getattr(n, "demand", None)
            if d is not None:
                T = int(np.asarray(d).reshape(-1).shape[0])
                break
        if T is None:
            for l in lines:
                m = getattr(l, "m", None)
                if m is not None and np.asarray(m).ndim >= 1:
                    T = int(np.asarray(m).reshape(-1).shape[0])
                    break
        if T is None:
            T = 1
    else:
        T = int(T)

    system.pressure_cache =  {
        "datum_node_id": int(datum_id),
        "datum_pos": int(datum_pos),
        "node_ids": node_ids,
        "id_to_pos": id_to_pos,
        "parent": parent,
        "parent_edge": parent_edge,
        "children": children,
        "order": np.array(order, dtype=int),
        "sign": sign,
        "L": L,
        "D": D,
        "eps": eps,
        "T": T,
    }              
                  
def _val_at_t(x, t, default=None):
    if x is None:
        return float(default)
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        return float(a)
    a = a.reshape(-1)
    return float(a[t]) if t < a.shape[0] else float(a[-1])


def compute_pressures_fast_tree(
    system,
    t,
    rho=1000.0,
    mu_supply=0.0003372,
    mu_return=0.0008891,
    return_pressure_default=150_000.0,
    required_pressure_difference_default=100_000.0,
    attr_return_p_set="return_pressure_setpoint",
    attr_required_dp="required_pressure_difference",
):
    nodes = system.nodes
    lines = system.lines
    n_nodes = len(nodes)
    n_lines = len(lines)
    cache = system.pressure_cache
    T = int(cache["T"])

    for n in nodes:
        if getattr(n, "supply_pressure", None) is None:
            n.supply_pressure = np.full(T, np.nan, dtype=float)
        if getattr(n, "return_pressure", None) is None:
            n.return_pressure = np.full(T, np.nan, dtype=float)

    for l in lines:
        if getattr(l, "supply_pressure_drop", None) is None:
            l.supply_pressure_drop = np.full(T, np.nan, dtype=float)
        if getattr(l, "return_pressure_drop", None) is None:
            l.return_pressure_drop = np.full(T, np.nan, dtype=float)

    m_vec = np.zeros(n_lines, dtype=float)
    for i, l in enumerate(lines):
        m_arr = getattr(l, "m_init", None)
        if m_arr is None:
            m_arr = getattr(l, "m", None)
        m_vec[i] = _val_at_t(m_arr, t, None)

    L = cache["L"]
    D = cache["D"]
    eps = cache["eps"] / 1000  
    A = np.pi * D**2 / 4.0

    v = m_vec / (rho * A)
    v_abs = np.abs(v)
    Re = np.maximum(rho * v_abs * D / mu_supply, 1e-12)
    f_lam = 64.0 / Re
    f_turb = 0.25 / (np.log10(eps/(3.7*D) + 5.74/(Re**0.9))**2)
    f = np.where(Re < 2300.0, f_lam, f_turb)
    dp_supply = f * (L/D) * 0.5 * rho * v * v_abs

    v_r = m_vec / (rho * A)
    v_r_abs = np.abs(v_r)
    Re_r = np.maximum(rho * v_r_abs * D / mu_return, 1e-12)
    f_lam_r = 64.0 / Re_r
    f_turb_r = 0.25 / (np.log10(eps/(3.7*D) + 5.74/(Re_r**0.9))**2)
    f_r = np.where(Re_r < 2300.0, f_lam_r, f_turb_r)
    dp_return = f_r * (L/D) * 0.5 * rho * v_r * v_r_abs

    for i, l in enumerate(lines):
        l.supply_pressure_drop[t] = float(dp_supply[i])
        l.return_pressure_drop[t] = float(dp_return[i])


    datum = int(cache["datum_pos"])
    parent = cache["parent"]
    parent_edge = cache["parent_edge"]
    sign = cache["sign"]
    order = cache["order"]


    Psup = np.full(n_nodes, np.nan, dtype=float)
    Psup[datum] = 0.0  

    for vpos in order:
        if vpos == datum:
            continue
        ppos = int(parent[vpos])
        eidx = int(parent_edge[vpos])
        Psup[vpos] = Psup[ppos] + sign[vpos] * dp_supply[eidx]

 
    Pret = np.full(n_nodes, np.nan, dtype=float)

    P_return0 = return_pressure_default
    p_ret_set = getattr(nodes[datum], attr_return_p_set, None)
    if p_ret_set is not None:
        P_return0 = _val_at_t(p_ret_set, t, return_pressure_default)

    Pret[datum] = float(P_return0)

    for vpos in order:
        if vpos == datum:
            continue
        ppos = int(parent[vpos])
        eidx = int(parent_edge[vpos])
        Pret[vpos] = Pret[ppos] - sign[vpos] * dp_return[eidx]


    required_dp = required_pressure_difference_default
    req_arr = getattr(nodes[datum], attr_required_dp, None)
    if req_arr is not None:
        required_dp = _val_at_t(req_arr, t, required_pressure_difference_default)

    delta = Psup - Pret
    min_delta = float(np.nanmin(delta))

    shift = float(required_dp) - min_delta
    Psup += shift  
    for i, n in enumerate(nodes):
        n.supply_pressure[t] = float(Psup[i])
        n.return_pressure[t] = float(Pret[i])



def _infer_T_from_nodes(system):
    for n in system.nodes:
        d = getattr(n, "demand", None)
        if d is not None:
            return int(np.asarray(d, dtype=float).reshape(-1).shape[0])


import numpy as np
from collections import deque


def initialize_mass_flow_guesses_tree_routing_timeseries(
    system,
    Tin_default=60.0,
    Tout_default=40.0,
    eff_default=0.99,
    cp=4.17,
    min_velocity=0.15,
    rho=1000.0,
    attr_tin="requested_inlet_temperature",
    attr_tout="expected_outlet_temperature",
    attr_eff="node_dhn_efficiency",
    seed=None,
):

    nodes = system.nodes
    lines = system.lines
    n_nodes = len(nodes)
    n_lines = len(lines)

    T = _infer_T_from_nodes(system)

    node_ids = np.array([n.node_id for n in nodes], dtype=int)
    max_id = node_ids.max() + 1
    id_to_pos = np.full(max_id, -1, dtype=int)
    id_to_pos[node_ids] = np.arange(n_nodes)

    rng = np.random.default_rng(seed)

    supply_mask = np.array(
        [getattr(n, "node_type", None) == "supply" for n in nodes],
        dtype=bool
    )
    supply_ids = node_ids[supply_mask]

    if supply_ids.size:
        datum_id = int(getattr(system, "datum_node",
                               rng.choice(supply_ids)))
    else:
        datum_id = int(node_ids[0])

    datum = id_to_pos[datum_id]

    consumers = [n for n in nodes
                 if getattr(n, "node_type", None) == "consumer"]

    consumer_positions = np.array(
        [id_to_pos[n.node_id] for n in consumers],
        dtype=int
    )

    tree_idx = np.asarray(system.spanning_tree_index, dtype=int)

    line_start = np.array([l.start_node for l in lines], dtype=int)
    line_end   = np.array([l.end_node   for l in lines], dtype=int)

    adj = [[] for _ in range(n_nodes)]

    for eidx in tree_idx:
        u = id_to_pos[line_start[eidx]]
        v = id_to_pos[line_end[eidx]]
        adj[u].append((v, eidx))
        adj[v].append((u, eidx))

    parent = np.full(n_nodes, -1, dtype=int)
    parent_edge = np.full(n_nodes, -1, dtype=int)

    parent[datum] = datum
    q = deque([datum])

    while q:
        u = q.popleft()
        for v, eidx in adj[u]:
            if parent[v] != -1:
                continue
            parent[v] = u
            parent_edge[v] = eidx
            q.append(v)

    consumer_paths = []

    for pos in consumer_positions:
        path_edges = []
        x = pos
        while x != datum and parent[x] != -1:
            path_edges.append(parent_edge[x])
            x = parent[x]
        consumer_paths.append(path_edges)

    D = np.array(
        [getattr(l, "inner_diameter", 0.0) for l in lines],
        dtype=float
    ) / 1000.0

    A = np.pi * D * D / 4.0
    m_min_line = rho * A * min_velocity
    m_min_line = np.where(np.isfinite(m_min_line) &
                          (m_min_line > 0.0),
                          m_min_line, 0.0)

    m_branch_all = np.zeros((T, n_lines), dtype=float)

    for t in range(T):

        dT = Tin_default - Tout_default
        if dT <= 0:
            continue

        denom = cp * dT * eff_default

        m_extract = np.zeros(n_nodes)

        for idx, n in enumerate(consumers):
            P = np.asarray(n.demand).reshape(-1)[t]
            m = P / denom if denom != 0 else 0.0
            if m < 0 or not np.isfinite(m):
                m = 0.0
            m_extract[consumer_positions[idx]] = m

        total_extract = m_extract.sum()

        for idx, m in enumerate(m_extract[consumer_positions]):
            if m == 0.0:
                continue
            for eidx in consumer_paths[idx]:
                m_branch_all[t, eidx] += m

        sgn = np.sign(m_branch_all[t])
        sgn[sgn == 0.0] = 1.0
        m_branch_all[t] = sgn * np.maximum(
            np.abs(m_branch_all[t]),
            m_min_line
        )

    for j, l in enumerate(lines):
        l.m_init = m_branch_all[:, j]

    for i, n in enumerate(nodes):
        n.m_init = np.zeros(T)

    system.init_guess_meta = {
        "method": "tree_routing_timeseries_fast",
        "datum_node_id": int(datum_id),
        "T": int(T),
    }

def earth_temperature_from_epw(
    epw_path: str,

    depth_m: float,
    k_soil: float = 1.5,
    rho_cp: float = 2.0e6,
) -> np.ndarray:
    """
    Returns hourly earth/soil temperature at a given depth (m) for an EPW year.

    Inputs
    - epw_path: path to EPW file
    - k_soil: soil thermal conductivity [W/m-K]
    - depth_m: depth below ground surface [m]
    - rho_cp: volumetric heat capacity of soil [J/m^3-K] (default ~ typical moist soil)

    Output
    - T_earth: np.ndarray shape (8760,), °C

    Model
    - Uses a 1D periodic conduction model forced by the *annual* sinusoid of air temperature:
        T(z,t) = T_mean + A*exp(-z/δ)*cos(ω*(t - t0) - z/δ)
      where δ = sqrt(2*α/ω), α = k/(rho_cp), ω = 2π / (365 days)
    - T_mean and A are estimated from the EPW dry-bulb.
    - t0 is aligned to the hour of minimum daily-mean air temperature.
    """
    if depth_m < 0:
        depth_m = abs(depth_m)

    alpha = float(k_soil) / float(rho_cp)  # m^2/s
    omega = 2.0 * np.pi / (365.0 * 24.0 * 3600.0)  # rad/s
    delta = np.sqrt(2.0 * alpha / omega) if alpha > 0 else np.inf

    with open(epw_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    data_lines = lines[8:] if len(lines) > 8 else []
    n = len(data_lines)

    T_air = np.full(8760, np.nan, dtype=float)

    for i in range(min(8760, n)):
        parts = data_lines[i].strip().split(",")
        if len(parts) < 7:
            continue
        try:
            T_air[i] = float(parts[6])  # Dry Bulb Temperature (C)
        except Exception:
            T_air[i] = np.nan

    if not np.any(np.isfinite(T_air)):
        return np.full(8760, np.nan, dtype=float)

    if T_air.size < 8760:
        tmp = np.full(8760, np.nan, dtype=float)
        tmp[: T_air.size] = T_air
        T_air = tmp

    if np.any(~np.isfinite(T_air)):
        idx = np.arange(T_air.size)
        good = np.isfinite(T_air)
        if good.sum() >= 2:
            T_air[~good] = np.interp(idx[~good], idx[good], T_air[good])
        else:
            fill = float(np.nanmean(T_air[good])) if good.sum() else 0.0
            T_air[:] = fill

    T_mean = float(np.mean(T_air))
    A0 = 0.5 * (float(np.max(T_air)) - float(np.min(T_air)))

    daily = T_air.reshape(-1, 24).mean(axis=1) if T_air.size >= 24 else np.array([T_mean])
    day_min = int(np.argmin(daily))
    t0_hours = float(day_min * 24)
    t0 = t0_hours * 3600.0

    t_sec = np.arange(8760, dtype=float) * 3600.0

    atten = np.exp(-depth_m / delta) if np.isfinite(delta) and delta > 0 else 0.0
    phase = depth_m / delta if np.isfinite(delta) and delta > 0 else 0.0

    T_earth = T_mean + (A0 * atten) * np.cos(omega * (t_sec - t0) - phase)
    return T_earth   


def _select_dn_index(D_required):
    idx = np.searchsorted(DN_INNER_D, D_required, side="left")
    if idx >= len(DN_INNER_D):
        idx = len(DN_INNER_D) - 1
    return idx

def rebuild_spanning_tree_index(system):
    from collections import deque

    nodes = system.nodes
    lines = system.lines

    node_ids = [n.node_id for n in nodes]
    id_to_pos = {nid: i for i, nid in enumerate(node_ids)}

    # find supply
    supply_ids = [n.node_id for n in nodes if n.node_type == "supply"]
    if not supply_ids:
        raise ValueError("No supply node")

    root = id_to_pos[supply_ids[0]]

    # adjacency with line indices
    adj = [[] for _ in nodes]
    for i, l in enumerate(lines):
        if not l.active:
            continue
        u = id_to_pos[l.start_node]
        v = id_to_pos[l.end_node]
        adj[u].append((v, i))
        adj[v].append((u, i))

    visited = [False] * len(nodes)
    visited[root] = True

    q = deque([root])
    spanning_tree_index = []

    while q:
        u = q.popleft()
        for v, eidx in adj[u]:
            if visited[v]:
                continue
            visited[v] = True
            spanning_tree_index.append(eidx)
            q.append(v)

    system.spanning_tree_index = spanning_tree_index


def lines_sizing(
    system,
    pipe_material="carbon_steel",
    insulation_class="series_2",
    delta_T=20.0,
    efficiency=0.99,
    delta_p_max=250.0,
    rho=1000.0,
    cp=4180.0,
    friction_factor=0.04,
    minimum_diameter_mm=5.0,
    concurrency="approx",
):

    material = PIPE_MATERIALS[pipe_material]
    insulation = INSULATION_CLASSES[insulation_class]


    for node in system.nodes:
        if getattr(node, "node_type", None) == "consumer":
            node.peak_demand = float(np.max(node.demand))
        else:
            node.peak_demand = 0.0

    node_ids = [n.node_id for n in system.nodes]
    id_to_pos = {nid: i for i, nid in enumerate(node_ids)}

    supply_nodes = [n.node_id for n in system.nodes if n.node_type == "supply"]
    if not supply_nodes:
        raise ValueError("No supply node found")

    root = id_to_pos[supply_nodes[0]]


    adj = [[] for _ in system.nodes]
    for eidx in system.spanning_tree_index:
        l = system.lines[eidx]
        u = id_to_pos[l.start_node]
        v = id_to_pos[l.end_node]
        adj[u].append((v, eidx))
        adj[v].append((u, eidx))

    parent = [-1] * len(system.nodes)
    parent_line = [-1] * len(system.nodes)
    parent[root] = root

    q = deque([root])
    order = []

    while q:
        u = q.popleft()
        order.append(u)
        for v, eidx in adj[u]:
            if parent[v] != -1:
                continue
            parent[v] = u
            parent_line[v] = eidx
            q.append(v)


    n_cons = np.zeros(len(system.nodes), dtype=int)
    raw_peak = np.zeros(len(system.nodes), dtype=float)

    for i, n in enumerate(system.nodes):
        if n.node_type == "consumer":
            n_cons[i] = 1
            raw_peak[i] = n.peak_demand

    for u in reversed(order):
        if u == root:
            continue
        p = parent[u]
        n_cons[p] += n_cons[u]
        raw_peak[p] += raw_peak[u]


    def g(n):
        if isinstance(concurrency, (int, float)):
            return float(concurrency)
        n = max(int(n), 1)
        return (65.0 * (1.2 ** (-0.05 * n)) + 35.0) / 100.0


    for l in system.lines:
        l.n_consumers = 0
        l.raw_peak = 0.0
        l.concurrent_peak = 0.0
        l.m_sizing = 0.0


    for v in range(len(system.nodes)):
        eidx = parent_line[v]
        if eidx < 0:
            continue

        nc = int(n_cons[v])
        rp = float(raw_peak[v])
        cp_peak = g(nc) * rp

        if cp_peak > 0:
            m = (cp_peak * 1000.0) / (cp * delta_T * efficiency)
        else:
            m = 0.0

        line = system.lines[eidx]
        L = float(line.length)

        if m > 0.0 and L > 0.0:
            D_req = (
                (8.0 * friction_factor * m * m)
                / (np.pi ** 2 * rho * delta_p_max)
            ) ** 0.2
        else:
            D_req = 0.0

        D_req = max(D_req, minimum_diameter_mm * 1e-3)

        dn_idx = _select_dn_index(D_req)
        line.inner_diameter_required = 0

        line.n_consumers = nc
        line.raw_peak = rp
        line.concurrent_peak = cp_peak
        line.m_sizing = m
        line.dn = DN_CATALOG[dn_idx]["dn"]
        line.inner_diameter = DN_INNER_D[dn_idx]*1000
        line.outer_diameter = DN_OUTER_D[dn_idx]*1000
        line.pipe_material = pipe_material
        line.pipe_conductivity = material["pipe_conductivity"]
        line.roughness = material["roughness"]
        line.pipe_specific_heat_capacity = material["pipe_specific_heat_capacity"]
        line.insulation_class = insulation_class
        line.insulation_thickness = insulation["insulation_thickness"]*1000
        line.insulation_conductivity = insulation["insulation_conductivity"]
        line.inner_diameter_required = D_req
        line.dn_index = dn_idx


def precompute_line_UA(
    system,
    burial_depth=1.0,
    soil_k_default=1.5,
    min_k=1e-6,
):
    for line in system.lines:
        D_i = line.inner_diameter / 1000
        D_o = line.outer_diameter / 1000

        t_ins = float(getattr(line, "insulation_thickness", 0.0))
        t_ins = t_ins / 1000.0 

        k_ins = float(getattr(line, "insulation_conductivity", 0.03))
        k_soil = soil_k_default
        k_ins = max(k_ins, min_k)
        k_soil = max(k_soil, min_k)
        # print(D_i)
        r_i = D_i / 2.0
        r_o = r_i + t_ins if t_ins > 0 else D_o / 2.0
        r_o = max(r_o, r_i * 1.0001)
        R_ins_per_m = np.log(r_o / r_i) / (2.0 * np.pi * k_ins)

        D_bury = 2.0 * r_o
        R_soil_per_m = np.log((4.0 * burial_depth) / max(D_bury, 1e-9)) / (2.0 * np.pi * k_soil)
        R_soil_per_m = max(R_soil_per_m, 0.0)

        R_tot_per_m = R_ins_per_m + R_soil_per_m
        UA = float(line.length) / max(R_tot_per_m, 1e-9)

        line.UA = UA
        

def _ensure_arr(obj, name, T):
    a = getattr(obj, name, None)
    if a is None:
        setattr(obj, name, np.full(T, np.nan, dtype=float))
        return
    a = np.asarray(a)
    if a.ndim == 0:
        setattr(obj, name, np.full(T, float(a), dtype=float))
        return


def solve_thermal_matrix_meanT_tree(
    system,
    t,
    init_supply_station_temperature=80.0,
    required_supply_temperature=60.0,
    threshold_temperature=2.0,
    cp=4180.0,
    eff_default=0.99,
    attr_eff="node_dhn_efficiency",
    demand_in_kw=True,
    eps_m=1e-9
):
    earth_temperature = system.ground_temperature_array[t]
    system.pressure_cache["earth_temperature"] = earth_temperature
    cache = system.pressure_cache
    nodes = system.nodes
    lines = system.lines
    id_to_pos = cache["id_to_pos"]
    parent = cache["parent"]
    parent_edge = cache["parent_edge"]
    children = cache["children"]
    Tlen = int(cache["T"])
    N = len(nodes)

    for n in nodes:
        _ensure_arr(n, "supply_temperature", Tlen)
        _ensure_arr(n, "return_temperature", Tlen)

    for l in lines:
        _ensure_arr(l, "supply_start_temperature", Tlen)
        _ensure_arr(l, "supply_end_temperature", Tlen)
        _ensure_arr(l, "return_start_temperature", Tlen)
        _ensure_arr(l, "return_end_temperature", Tlen)

    supply_pos = [i for i, n in enumerate(nodes) if getattr(n, "node_type", None) == "supply"]
    cons_pos = [i for i, n in enumerate(nodes) if getattr(n, "node_type", None) == "consumer"]

    def _solve_once(T_supply):
        m_edge = np.zeros(len(lines), dtype=float)
        a_edge = np.zeros(len(lines), dtype=float)
        b_edge = np.zeros(len(lines), dtype=float)

        for eidx, l in enumerate(lines):
            m_arr = getattr(l, "m_init", None)
            if m_arr is None:
                m_arr = getattr(l, "m", None)
            if m_arr is None:
                m = 0.0
            else:
                m = _val_at_t(m_arr, t, 0.0)

            mdot = max(abs(float(m)), float(eps_m))
            U = float(getattr(l, "UA", 0.0))
            Mflow = mdot * float(cp)
            denom = (Mflow + 0.5 * U)
            if abs(denom) < 1e-30:
                a = 1.0
                b = 0.0
            else:
                a = (Mflow - 0.5 * U) / denom
                b = (U * float(earth_temperature)) / denom

            m_edge[eidx] = mdot
            a_edge[eidx] = a
            b_edge[eidx] = b

        def idxS(i): return i
        def idxR(i): return N + i

        rows = []
        cols = []
        data = []
        rhs = np.zeros(2 * N, dtype=float)
        r = 0

        supply_set = set(supply_pos)

        for i in range(N):
            if i in supply_set:
                rows.append(r); cols.append(idxS(i)); data.append(1.0)
                rhs[r] = float(T_supply)
                r += 1
            else:
                p = int(parent[i])
                eidx = int(parent_edge[i])
                a = float(a_edge[eidx])
                b = float(b_edge[eidx])

                rows += [r, r]
                cols += [idxS(i), idxS(p)]
                data += [1.0, -a]
                rhs[r] = b
                r += 1

        cons_set = set(cons_pos)

        for i in range(N):
            n = nodes[i]
            if i in cons_set:
                P = _val_at_t(getattr(n, "demand", None), t, 0.0)
                if demand_in_kw:
                    P = float(P) * 1000.0
                eff = _val_at_t(getattr(n, attr_eff, None), t, eff_default)
                eff = max(float(eff), 1e-9)

                mnode_arr = getattr(n, "m_init", None)
                mdot_node = max(abs(_val_at_t(mnode_arr, t, 0.0)), float(eps_m))

                dT = float(P) / (eff * mdot_node * float(cp)) if mdot_node > 0 else 0.0
                rows += [r, r]
                cols += [idxR(i), idxS(i)]
                data += [1.0, -1.0]
                rhs[r] = -float(dT)
                r += 1
            else:
                ch = children[i]
                if len(ch) == 0:
                    rows += [r, r]
                    cols += [idxR(i), idxS(i)]
                    data += [1.0, -1.0]
                    rhs[r] = 0.0
                    r += 1
                else:
                    W = 0.0
                    for c in ch:
                        eidx = int(parent_edge[c])
                        W += float(m_edge[eidx])
                    W = max(W, float(eps_m))

                    rows.append(r); cols.append(idxR(i)); data.append(1.0)

                    rhs_mix = 0.0
                    for c in ch:
                        eidx = int(parent_edge[c])
                        w = float(m_edge[eidx])
                        a = float(a_edge[eidx])
                        b = float(b_edge[eidx])
                        rows.append(r); cols.append(idxR(int(c))); data.append(-(w * a) / W)
                        rhs_mix += (w * b) / W

                    rhs[r] = float(rhs_mix)
                    r += 1

        if _HAS_SCIPY:
            try:
                A = coo_matrix((data, (rows, cols)), shape=(2 * N, 2 * N)).tocsr()
                x = spsolve(A, rhs)
            except Exception:
                x = np.full(2 * N, np.nan, dtype=float)
        else:
            try:
                A = np.zeros((2 * N, 2 * N), dtype=float)
                for rr, cc, vv in zip(rows, cols, data):
                    A[int(rr), int(cc)] += float(vv)
                x = np.linalg.solve(A, rhs)
            except Exception:
                x = np.full(2 * N, np.nan, dtype=float)

        Ts = np.asarray(x[:N], dtype=float)
        Tr = np.asarray(x[N:], dtype=float)
        return Ts, Tr

    def _pick_first_supply_guess():
        if not supply_pos:
            return float(init_supply_station_temperature)

        p = int(supply_pos[0])
        try:
            arr = getattr(nodes[p], "supply_temperature", None)
            if arr is None:
                return float(init_supply_station_temperature)
            v = float(np.asarray(arr, dtype=float).reshape(-1)[t])
            if np.isfinite(v):
                return v
        except Exception:
            pass
        return float(init_supply_station_temperature)

    def _worst_Ts(Ts):
        if cons_pos:
            vals = np.asarray([Ts[i] for i in cons_pos], dtype=float)
            if np.any(np.isfinite(vals)):
                j = int(np.nanargmin(vals))
                i = int(cons_pos[j])
                return float(Ts[i]), i
        if np.any(np.isfinite(Ts)):
            i = int(np.nanargmin(np.asarray(Ts, dtype=float)))
            return float(Ts[i]), i
        return float("nan"), 0
    
    
    ok_band = False
    if hasattr(system, "last_supply_temperature"):
        T_supply = system.last_supply_temperature 
    else:
        T_supply = init_supply_station_temperature
        
    lo = float(required_supply_temperature)
    hi = float(required_supply_temperature) + float(threshold_temperature)   
    while not(ok_band):
        Ts, Tr = _solve_once(T_supply)
        worst, worst_pos = _worst_Ts(Ts)
        ok_band = (worst >= lo) and (worst <= hi)
        if ok_band:
            shift = 0.0 
            T_supply = T_supply + shift 
        else:
            if worst < lo:
                shift = hi - worst 
            if worst > hi :
                shift = lo - worst 
            
            T_supply = T_supply + shift 
        
        system.last_supply_temperature = T_supply 

        
        
    for i, n in enumerate(nodes):
        n.supply_temperature[t] = float(Ts[i])
        n.return_temperature[t] = float(Tr[i])

    for l in lines:
        try:
            s_id = int(l.start_node)
            e_id = int(l.end_node)
        except Exception:
            continue
        if s_id not in id_to_pos or e_id not in id_to_pos:
            continue
        s = int(id_to_pos[s_id])
        e = int(id_to_pos[e_id])
        l.supply_start_temperature[t] = float(Ts[s])
        l.supply_end_temperature[t] = float(Ts[e])
        l.return_start_temperature[t] = float(Tr[s])
        l.return_end_temperature[t] = float(Tr[e])

    cache["node_with_lowest_supply_temperature"] = int(getattr(nodes[int(worst_pos)], "node_id", int(worst_pos) + 1))
    cache["node_with_lowest_supply_temperature_pos"] = int(worst_pos)

    system.meta_last_thermal_step = {
        "timestep": int(t),
        "method": "matrix_meanT_tree_shifted",
        "T_ground": float(earth_temperature),
        "cp": float(cp),
        "eff_default": float(eff_default),
        "T_supply_first": float(T_supply),
        "worst_supply_temp": float(worst) if np.isfinite(worst) else None,
        "ok_in_band": bool(ok_band),
        "node_with_lowest_supply_temperature": int(cache["node_with_lowest_supply_temperature"]),
    }

    return Ts, Tr



def add_line_thermal_losses(system, cp=4180.0, eps_m=1e-9):
    def _ensure_arr(obj, name: str, T: int):
        a = getattr(obj, name, None)
        if a is None:
            setattr(obj, name, np.full(T, np.nan, dtype=float))
            return
        a = np.asarray(a, dtype=float)
        if a.ndim == 0:
            setattr(obj, name, np.full(T, float(a), dtype=float))
            return
        a = a.reshape(-1)
        if a.shape[0] != T:
            b = np.full(T, np.nan, dtype=float)
            n = min(T, a.shape[0])
            b[:n] = a[:n]
            setattr(obj, name, b)

    T = None
    for l in system.lines:
        a = getattr(l, "supply_start_temperature", None)
        if a is not None:
            T = int(np.asarray(a).reshape(-1).shape[0])
            break
    if T is None:
        raise ValueError("Cannot infer T. Run thermal solve first so line temperature arrays exist.")

    cpj = float(cp) if float(cp) > 100.0 else float(cp) * 1000.0

    for l in system.lines:
        _ensure_arr(l, "supply_loss", T)
        _ensure_arr(l, "return_loss", T)

        Tin_s = np.asarray(l.supply_start_temperature, dtype=float).reshape(-1)
        Tout_s = np.asarray(l.supply_end_temperature, dtype=float).reshape(-1)
        Tin_r = np.asarray(l.return_start_temperature, dtype=float).reshape(-1)
        Tout_r = np.asarray(l.return_end_temperature, dtype=float).reshape(-1)

        m = getattr(l, "m_init", None)
        if m is None:
            m = getattr(l, "m", None)
        if m is None:
            raise ValueError(f"Line {getattr(l,'line_id','?')} missing 'm_init' or 'm'.")
        m = np.asarray(m, dtype=float).reshape(-1)
        if m.shape[0] == 1:
            mdot = np.full(T, float(m[0]), dtype=float)
        else:
            mdot = m[:T]

        mdot_abs = np.maximum(np.abs(mdot), eps_m)

        l.supply_loss[:T] = mdot_abs * cpj * np.maximum(Tin_s[:T] - Tout_s[:T], 0.0)
        l.return_loss[:T] = mdot_abs * cpj * np.maximum(Tin_r[:T] - Tout_r[:T], 0.0)

def solve_thermal_matrix_meanT_tree_no_loss(
    system,
    t,
    required_supply_temperature= 60,
    cp=4180.0,
    eff_default=0.99,
    attr_eff="node_dhn_efficiency",
    demand_in_kw=True,
    eps_m=1e-9,
):
    """
    Tree thermal solve with NO distribution losses (UA=0).
    Writes results into *_no_loss attributes:

    Nodes:
      - supply_temperature_no_loss
      - return_temperature_no_loss

    Lines:
      - supply_start_temperature_no_loss
      - supply_end_temperature_no_loss
      - return_start_temperature_no_loss
      - return_end_temperature_no_loss
    """
    cache = system.pressure_cache
    nodes = system.nodes
    lines = system.lines
    id_to_pos = cache["id_to_pos"]
    parent = cache["parent"]
    parent_edge = cache["parent_edge"]
    children = cache["children"]
    Tlen = int(cache["T"])
    N = len(nodes)


    for n in nodes:
        _ensure_arr(n, "supply_temperature_no_loss", Tlen)
        _ensure_arr(n, "return_temperature_no_loss", Tlen)

    for l in lines:
        _ensure_arr(l, "supply_start_temperature_no_loss", Tlen)
        _ensure_arr(l, "supply_end_temperature_no_loss", Tlen)
        _ensure_arr(l, "return_start_temperature_no_loss", Tlen)
        _ensure_arr(l, "return_end_temperature_no_loss", Tlen)

    supply_pos = [i for i, n in enumerate(nodes) if getattr(n, "node_type", None) == "supply"]
    cons_pos = [i for i, n in enumerate(nodes) if getattr(n, "node_type", None) == "consumer"]

    supply_set = set(supply_pos)
    cons_set = set(cons_pos)

    m_edge = np.zeros(len(lines), dtype=float)
    for eidx, l in enumerate(lines):
        m_arr = getattr(l, "m", None)
        if m_arr is None:
            m_arr = getattr(l, "m_init", None)
        if m_arr is None:
            m_edge[eidx] = eps_m
        else:
            m_edge[eidx] = max(abs(_val_at_t(m_arr, t, 0.0)), eps_m)

    def idxS(i): return i
    def idxR(i): return N + i

    rows, cols, data = [], [], []
    rhs = np.zeros(2 * N, dtype=float)
    r = 0

    for i in range(N):
        if i in supply_set:
            rows.append(r); cols.append(idxS(i)); data.append(1.0)
            rhs[r] = float(required_supply_temperature)
            r += 1
        else:
            p = int(parent[i])
            rows += [r, r]
            cols += [idxS(i), idxS(p)]
            data += [1.0, -1.0]
            rhs[r] = 0.0
            r += 1

    for i in range(N):
        n = nodes[i]
        if i in cons_set:
            P = _val_at_t(getattr(n, "demand", None), t, 0.0)
            if demand_in_kw:
                P *= 1000.0

            eff = _val_at_t(getattr(n, attr_eff, None), t, eff_default)
            eff = max(float(eff), 1e-9)

            mnode_arr = getattr(n, "m_init", None)
            mdot_node = max(abs(_val_at_t(mnode_arr, t, 0.0)), eps_m)

            dT = float(P) / (eff * mdot_node * float(cp))*0.5

            rows += [r, r]
            cols += [idxR(i), idxS(i)]
            data += [1.0, -1.0]
            rhs[r] = -dT
            r += 1
        else:
            ch = children[i]
            if len(ch) == 0:
                rows += [r, r]
                cols += [idxR(i), idxS(i)]
                data += [1.0, -1.0]
                rhs[r] = 0.0
                r += 1
            else:
                W = 0.0
                weights = []
                for c in ch:
                    eidx = int(parent_edge[c])
                    w = float(m_edge[eidx])
                    weights.append((int(c), w))
                    W += w
                W = max(W, eps_m)

                rows.append(r); cols.append(idxR(i)); data.append(1.0)
                rhs[r] = 0.0

                for c, w in weights:
                    rows.append(r); cols.append(idxR(c)); data.append(-(w / W))

                r += 1

    if _HAS_SCIPY:
        A = coo_matrix((data, (rows, cols)), shape=(2 * N, 2 * N)).tocsr()
        x = spsolve(A, rhs)
    else:
        A = np.zeros((2 * N, 2 * N), dtype=float)
        for rr, cc, vv in zip(rows, cols, data):
            A[rr, cc] += vv
        x = np.linalg.solve(A, rhs)

    Ts = x[:N].astype(float)
    Tr = x[N:].astype(float)

    for i, n in enumerate(nodes):
        n.supply_temperature_no_loss[t] = float(Ts[i])
        n.return_temperature_no_loss[t] = float(Tr[i])

    for l in lines:
        s_id = int(l.start_node)
        e_id = int(l.end_node)
        if s_id not in id_to_pos or e_id not in id_to_pos:
            continue
        s = id_to_pos[s_id]
        e = id_to_pos[e_id]
        l.supply_start_temperature_no_loss[t] = float(Ts[s])
        l.supply_end_temperature_no_loss[t] = float(Ts[e])
        l.return_start_temperature_no_loss[t] = float(Tr[s])
        l.return_end_temperature_no_loss[t] = float(Tr[e])

    return Ts, Tr

def _infer_T_from_node(system, attr):
    for n in system.nodes:
        a = getattr(n, attr, None)
        if a is not None:
            a = np.asarray(a).reshape(-1)
            if a.size > 1:
                return int(a.size)
            return 1
def _build_distance_from_source(system, cache):
    node_ids = [int(n.node_id) for n in system.nodes]
    pos_by_id = {nid: i for i, nid in enumerate(node_ids)}

    parent = np.asarray(cache["parent"], dtype=int)
    parent_edge = np.asarray(cache["parent_edge"], dtype=int)
    root = int(cache["datum_pos"])

    L = np.array([float(l.length) for l in system.lines], dtype=float)

    dist = np.full(len(node_ids), np.nan, dtype=float)
    dist[root] = 0.0

    order = cache.get("order", None)
    if order is None:
        order = list(range(len(node_ids)))

    changed = True
    it = 0
    while changed and it < len(node_ids) + 5:
        changed = False
        it += 1
        for i in order:
            if i == root:
                continue
            p = parent[i]
            eidx = parent_edge[i]
            if p < 0 or eidx < 0:
                continue
            if np.isfinite(dist[p]) and not np.isfinite(dist[i]):
                dist[i] = dist[p] + L[eidx]
                changed = True

    return dist, pos_by_id, root


def _node_vec(system, attr, t):
    v = np.array([float(np.asarray(getattr(n, attr)).reshape(-1)[t]) for n in system.nodes], dtype=float)
    return v

def add_hourly_energy_and_pump_metrics(system, rho=1000.0, cp=4180.0, pump_eff=0.7):
    def _as1d(x):
        a = np.asarray(x, dtype=float)
        if a.ndim == 0:
            return a.reshape(1)
        return a.reshape(-1)

    supply_nodes = [n for n in system.nodes if getattr(n, "node_type", None) == "supply"]
    if not supply_nodes:
        raise ValueError("No supply node in system.")
    supply = supply_nodes[0]
    sid = int(supply.node_id)

    Ts = getattr(supply, "supply_temperature", None)
    Tr = getattr(supply, "return_temperature", None)
    Ps = getattr(supply, "supply_pressure", None)
    Pr = getattr(supply, "return_pressure", None)

    if Ts is None or Tr is None:
        raise ValueError("Supply node missing supply_temperature/return_temperature.")
    if Ps is None or Pr is None:
        raise ValueError("Supply node missing supply_pressure/return_pressure.")

    Ts = _as1d(Ts)
    Tr = _as1d(Tr)
    Ps = _as1d(Ps)
    Pr = _as1d(Pr)

    T = int(min(Ts.shape[0], Tr.shape[0], Ps.shape[0], Pr.shape[0]))

    cpj = float(cp) if float(cp) > 100.0 else float(cp) * 1000.0

    hourly_heat_demand = np.zeros(T, dtype=float)
    for n in system.nodes:
        if getattr(n, "node_type", None) == "consumer":
            d = getattr(n, "demand", None)
            if d is None:
                continue
            d = _as1d(d)
            if d.shape[0] == 1:
                hourly_heat_demand += float(d[0]) * 1000.0
            else:
                hourly_heat_demand += d[:T] * 1000.0

    hourly_heat_loss = np.zeros(T, dtype=float)
    for l in system.lines:
        sl = getattr(l, "supply_loss", None)
        rl = getattr(l, "return_loss", None)
        if sl is not None:
            hourly_heat_loss += _as1d(sl)[:T]
        if rl is not None:
            hourly_heat_loss += _as1d(rl)[:T]

    mdot_out = np.zeros(T, dtype=float)
    for l in system.lines:
        m = getattr(l, "m_init", None)
        if m is None:

            m = getattr(l, "m", None)
        if m is None:
            continue
        m = _as1d(m)
        if m.shape[0] == 1:
            mt = np.full(T, float(m[0]), dtype=float)
        else:
            mt = m[:T]
        
        if int(l.start_node) == sid:
            flow_out = mt
        elif int(l.end_node) == sid:
            flow_out = -mt
        else:
            continue
        mdot_out += np.maximum(-flow_out, 0.0)

    heat_injected = mdot_out * cpj * np.maximum(Ts[:T] - Tr[:T], 0.0)

    dp = np.maximum(Ps[:T] - Pr[:T], 0.0)
    Vdot = mdot_out / float(rho)
    pump_power = dp * Vdot / float(pump_eff)

    system.hourly_heat_demand = hourly_heat_demand
    system.hourly_heat_loss = hourly_heat_loss

    supply.heat_injected = heat_injected
    supply.pump_power = pump_power

    return hourly_heat_demand, hourly_heat_loss, heat_injected, pump_power


def plot_info(
    system,
    name,
    cache=None,
    t0=1,
    start_datetime="2025-01-01 00:00:00",
    rho=1000.0,
    g=9.81,
    cp=4180.0,
    demand_in_kw=True,
    temp_with_loss_sup="supply_temperature",
    temp_with_loss_ret="return_temperature",
    temp_no_loss_sup="supply_temperature_no_loss",
    temp_no_loss_ret="return_temperature_no_loss",
    pressure_sup="supply_pressure",
    pressure_ret="return_pressure",
    mass_attr="m_init",
    mass_fallback="m",
    loss_attr_supply="supply_loss",
    loss_attr_return="return_loss",
    demand_attr="demand",
    consumer_types=("consumer",),
    color_supply="red",
    color_return="blue",
    color_loss="#7EC8E3",
    color_demand="#F4A261",
):



    if cache is None:
        cache = system.pressure_cache

    dist, pos_by_id, root_pos = _build_distance_from_source(system, cache)
    Lmax = float(np.nanmax(dist))

    T1 = _infer_T_from_node(system, temp_with_loss_sup)
    T2 = _infer_T_from_node(system, temp_no_loss_sup)
    T3 = _infer_T_from_node(system, pressure_sup)
    T_use = int(min(T1, T2, T3, int(cache.get("T", min(T1, T2, T3)))))

    seg_s = []
    seg_e = []
    plot_lines = []
    seg_a = []
    seg_b = []
    seg_len = []
    for l in system.lines:
        s_id = int(l.start_node)
        e_id = int(l.end_node)
        if s_id not in pos_by_id or e_id not in pos_by_id:
            continue
        s = pos_by_id[s_id]
        e = pos_by_id[e_id]
        seg_s.append(s)
        seg_e.append(e)
        plot_lines.append(l)

        a = float(min(dist[s], dist[e]))
        b = float(max(dist[s], dist[e]))
        L = max(b - a, 1e-12)
        seg_a.append(a)
        seg_b.append(b)
        seg_len.append(L)

    seg_s = np.asarray(seg_s, dtype=int)
    seg_e = np.asarray(seg_e, dtype=int)
    xs_sup = np.vstack([dist[seg_s], dist[seg_e]]).T
    xs_ret = np.vstack([dist[seg_e], dist[seg_s]]).T

    seg_a = np.asarray(seg_a, dtype=float)
    seg_b = np.asarray(seg_b, dtype=float)
    seg_len = np.asarray(seg_len, dtype=float)

    cons_pos = []
    cons_x = []
    for i, n in enumerate(system.nodes):
        if getattr(n, "node_type", None) in consumer_types:
            cons_pos.append(i)
            cons_x.append(float(dist[i]))
    cons_pos = np.asarray(cons_pos, dtype=int)
    cons_x = np.asarray(cons_x, dtype=float)

    def _line_m_at_t(line, t):
        arr = getattr(line, mass_attr, None)
        if arr is None:
            arr = getattr(line, mass_fallback, None)
        if arr is None:
            return 0.0
        a = np.asarray(arr, dtype=float)
        if a.ndim == 0:
            return float(a)
        a = a.reshape(-1)
        return float(a[t]) if t < a.shape[0] else float(a[-1])

    mmax_global = 1.0
    if plot_lines:
        vals = []
        for l in plot_lines:
            arr = getattr(l, mass_attr, None)
            if arr is None:
                arr = getattr(l, mass_fallback, None)
            if arr is None:
                continue
            a = np.asarray(arr, dtype=float)
            if a.ndim == 0:
                vals.append(abs(float(a)))
            else:
                a = a.reshape(-1)[:T_use]
                vals.append(float(np.nanmax(np.abs(a))) if a.size else 0.0)
        if vals:
            mmax_global = max(vals)
            if not np.isfinite(mmax_global) or mmax_global <= 0:
                mmax_global = 1.0
    mmax_global *= 1.1

    supply_losses = []
    return_losses = []
    for l in plot_lines:
        ls = getattr(l, loss_attr_supply, None)/1000
        lr = getattr(l, loss_attr_return, None)/1000
        ls = np.asarray(ls, dtype=float).reshape(-1)[:T_use]
        lr = np.asarray(lr, dtype=float).reshape(-1)[:T_use]
        supply_losses.append(ls)
        return_losses.append(lr)
    supply_losses = np.vstack(supply_losses) if supply_losses else np.zeros((0, T_use), dtype=float)
    return_losses = np.vstack(return_losses) if return_losses else np.zeros((0, T_use), dtype=float)

    demands = []
    for idx in cons_pos:
        d = getattr(system.nodes[int(idx)], demand_attr, None)
        if d is None:
            demands.append(np.zeros(T_use, dtype=float))
        else:
            a = np.asarray(d, dtype=float).reshape(-1)
            if a.size >= T_use:
                demands.append(a[:T_use])
            else:
                demands.append(np.pad(a, (0, T_use - a.size), mode="edge"))
    demands = np.vstack(demands) if demands else np.zeros((0, T_use), dtype=float)
    if demand_in_kw:
        demands = demands

    ex_out = np.concatenate([seg_a, seg_b]) if seg_a.size else np.zeros(0, dtype=float)
    ex_ret = np.concatenate([2.0 * Lmax - seg_b, 2.0 * Lmax - seg_a]) if seg_a.size else np.zeros(0, dtype=float)
    bp_static = np.array(sorted(set([0.0, Lmax, 2.0 * Lmax] + ex_out.tolist() + ex_ret.tolist() + cons_x.tolist())), dtype=float)
    bp_static = bp_static[(bp_static >= 0.0) & (bp_static <= 2.0 * Lmax)]
    if bp_static.size < 2:
        bp_static = np.array([0.0, 2.0 * Lmax], dtype=float)

    total_end = (np.nansum(supply_losses + return_losses, axis=0) + np.nansum(demands, axis=0))
    yheat_max = float(np.nanmax(total_end))
    if not np.isfinite(yheat_max) or yheat_max <= 0:
        yheat_max = 1.0
    yheat_max *= 1.1

    def _curve_at_t(t: int):
        t = int(t)
        slopes_out = (supply_losses[:, t] / seg_len) if supply_losses.size else np.zeros(seg_len.shape[0], dtype=float)
        slopes_ret = (return_losses[:, t] / seg_len) if return_losses.size else np.zeros(seg_len.shape[0], dtype=float)

        events_x = []
        events_ds = []

        if seg_a.size:
            for a, b, s in zip(seg_a, seg_b, slopes_out):
                if np.isfinite(s) and s != 0.0:
                    events_x.extend([float(a), float(b)])
                    events_ds.extend([+float(s), -float(s)])
            for a, b, s in zip(seg_a, seg_b, slopes_ret):
                if np.isfinite(s) and s != 0.0:
                    ra = float(2.0 * Lmax - b)
                    rb = float(2.0 * Lmax - a)
                    events_x.extend([ra, rb])
                    events_ds.extend([+float(s), -float(s)])

        dem_steps = {}
        if demands.size:
            dem_t = demands[:, t]
            for x, q in zip(cons_x, dem_t):
                if np.isfinite(q) and q != 0.0:
                    dem_steps[float(x)] = dem_steps.get(float(x), 0.0) + float(q)

        if events_x:
            ex = np.asarray(events_x, dtype=float)
            eds = np.asarray(events_ds, dtype=float)
            o = np.argsort(ex)
            ex = ex[o]
            eds = eds[o]
        else:
            ex = np.zeros(0, dtype=float)
            eds = np.zeros(0, dtype=float)

        y_loss = np.zeros(bp_static.shape[0], dtype=float)
        y_dem = np.zeros(bp_static.shape[0], dtype=float)

        slope = 0.0
        loss = 0.0
        demcum = 0.0
        ei = 0

        x_prev = float(bp_static[0])
        while ei < ex.size and ex[ei] == x_prev:
            slope += float(eds[ei])
            ei += 1
        if x_prev in dem_steps:
            demcum += float(dem_steps[x_prev])

        y_loss[0] = loss
        y_dem[0] = demcum

        for k in range(1, bp_static.shape[0]):
            x = float(bp_static[k])
            dx = x - x_prev
            if dx < 0:
                dx = 0.0
            loss += slope * dx
            while ei < ex.size and ex[ei] == x:
                slope += float(eds[ei])
                ei += 1
            if x in dem_steps:
                demcum += float(dem_steps[x])
            y_loss[k] = loss
            y_dem[k] = demcum
            x_prev = x

        return bp_static, y_loss, y_dem

    datum = int(cache.get("datum_pos", root_pos))
    if not (0 <= datum < len(system.nodes)):
        datum = root_pos

    supply_node = system.nodes[datum]

    def _safe_series(obj, name):
        a = getattr(obj, name, None)
        if a is None:
            return None
        a = np.asarray(a, dtype=float)
        if a.ndim == 0:
            return np.full(T_use, float(a), dtype=float)
        a = a.reshape(-1)
        if a.size < T_use:
            a = np.pad(a, (0, T_use - a.size), mode="edge")
        return a[:T_use]

    heat_inj = _safe_series(supply_node, "heat_injected")
    heat_loss = _safe_series(system, "hourly_heat_loss")
    pump_pow = _safe_series(supply_node, "pump_power")
    if heat_inj is None or pump_pow is None:
        ch = cache.get("children", None)
        pe = cache.get("parent_edge", None)
        first_edge_idx = None
        if ch is not None and pe is not None and ch[datum]:
            c0 = int(ch[datum][0])
            first_edge_idx = int(pe[c0])
        if first_edge_idx is None:
            first_edge_idx = 0
        l0 = system.lines[first_edge_idx]
        m0 = np.array([abs(_line_m_at_t(l0, t)) for t in range(T_use)], dtype=float)
        Ts0 = _safe_series(supply_node, temp_with_loss_sup)
        Tr0 = _safe_series(supply_node, temp_with_loss_ret)
        if Ts0 is None or Tr0 is None:
            Ts0 = np.zeros(T_use, dtype=float)
            Tr0 = np.zeros(T_use, dtype=float)
        heat_inj = m0 * float(cp) * (Ts0 - Tr0)

        Ps0 = _safe_series(supply_node, pressure_sup)
        Pr0 = _safe_series(supply_node, pressure_ret)
        if Ps0 is None or Pr0 is None:
            Ps0 = np.zeros(T_use, dtype=float)
            Pr0 = np.zeros(T_use, dtype=float)
        dp_pump = np.maximum(Ps0 - Pr0, 0.0)
        pump_pow = dp_pump * (m0 / float(rho))

    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(3, 5, figure=fig, wspace=0.35, hspace=0.5)

    axT = fig.add_subplot(gs[0, 0:2])
    axTnl = fig.add_subplot(gs[0, 2:4])
    axP = fig.add_subplot(gs[1, 0:2])
    axM = fig.add_subplot(gs[1, 2:4])
    axH = fig.add_subplot(gs[2, 0:4])

    axInfo1 = fig.add_subplot(gs[0, 4])
    axInfo2 = fig.add_subplot(gs[1, 4])
    axInfo3 = fig.add_subplot(gs[2, 4])

    for ax in (axInfo1, axInfo2, axInfo3):
        ax.set_axis_off()

    s_ax = fig.add_axes([0.965, 0.12, 0.015, 0.76])
    sld = Slider(s_ax, "t", 0, T_use - 1, valinit=int(t0), valstep=1, orientation="vertical")

    def _init_profile_ax(ax, title, ylab, ylim=None):
        ax.set_title(title)
        ax.set_xlabel("Length from source [m]")
        ax.set_ylabel(ylab)
        ax.grid(True)
        ax.set_xlim(0, float(np.nanmax(dist)*1.2))
        if ylim is not None:
            ax.set_ylim(*ylim)

    _init_profile_ax(axT, "Temperature (with loss)", "Temperature [°C]", (30, 100))
    _init_profile_ax(axTnl, "Temperature (no loss)", "Temperature [°C]", (30, 100))
    _init_profile_ax(axP, "Pressure", "Pressure [bar]", (1, 5))
    _init_profile_ax(axM, "Mass flow (lines)", "Mass flow [kg/s]", (-mmax_global, mmax_global))

    axH.set_title("Rejected Heat")
    axH.set_xlabel("Travelled distance")
    axH.set_ylabel("Rejected Heat [kWh]")
    axH.grid(True)
    axH.set_xlim(0.0, 2.0 * Lmax)
    axH.set_ylim(0.0, yheat_max)
    Tg_arr = np.asarray(system.ground_temperature_array, dtype=float).reshape(-1)
    
    Tg0 = Tg_arr[int(t0)] if Tg_arr.size else np.nan
    ground_line = axT.axhline(
        Tg0,
        color="brown",
        linestyle="--",
        linewidth=1.5,
        label="Ground temperature"
    )

    def _node_vec_attr(attr, t):
        return np.array([float(np.asarray(getattr(n, attr)).reshape(-1)[t]) for n in system.nodes], dtype=float)

    def _text(ax, s):
        ax.clear()
        ax.set_axis_off()
        ax.text(0.0, 1.0, s, va="top", ha="left", fontsize=10, family="monospace")

    t = int(t0)

    Ts = _node_vec_attr(temp_with_loss_sup, t)
    Tr = _node_vec_attr(temp_with_loss_ret, t)
    Ts_nl = _node_vec_attr(temp_no_loss_sup, t)
    Tr_nl = _node_vec_attr(temp_no_loss_ret, t)

    Ps = _node_vec_attr(pressure_sup, t) / 100_000.0
    Pr = _node_vec_attr(pressure_ret, t) / 100_000.0

    linesT_sup = []
    linesT_ret = []
    for i in range(xs_sup.shape[0]):
        lsup, = axT.plot(xs_sup[i], [Ts[seg_s[i]], Ts[seg_e[i]]], color=color_supply, linewidth=1.0)
        lret, = axT.plot(xs_ret[i], [Tr[seg_e[i]], Tr[seg_s[i]]], color=color_return, linewidth=1.0)
        linesT_sup.append(lsup)
        linesT_ret.append(lret)

    scT_sup = axT.scatter(dist, Ts, s=12, color=color_supply)
    scT_ret = axT.scatter(dist, Tr, s=12, color=color_return)

    linesTnl_sup = []
    linesTnl_ret = []
    for i in range(xs_sup.shape[0]):
        lsup, = axTnl.plot(xs_sup[i], [Ts_nl[seg_s[i]], Ts_nl[seg_e[i]]], color=color_supply, linewidth=1.0)
        lret, = axTnl.plot(xs_ret[i], [Tr_nl[seg_e[i]], Tr_nl[seg_s[i]]], color=color_return, linewidth=1.0)
        linesTnl_sup.append(lsup)
        linesTnl_ret.append(lret)

    scTnl_sup = axTnl.scatter(dist, Ts_nl, s=12, color=color_supply)
    scTnl_ret = axTnl.scatter(dist, Tr_nl, s=12, color=color_return)

    linesP_sup = []
    linesP_ret = []
    for i in range(xs_sup.shape[0]):
        lsup, = axP.plot(xs_sup[i], [Ps[seg_s[i]], Ps[seg_e[i]]], color=color_supply, linewidth=1.0)
        lret, = axP.plot(xs_ret[i], [Pr[seg_e[i]], Pr[seg_s[i]]], color=color_return, linewidth=1.0)
        linesP_sup.append(lsup)
        linesP_ret.append(lret)

    scP_sup = axP.scatter(dist, Ps, s=12, color=color_supply)
    scP_ret = axP.scatter(dist, Pr, s=12, color=color_return)

    m_abs_plot = np.asarray([abs(_line_m_at_t(l, t)) for l in plot_lines], dtype=float)
    linesM_sup = []
    linesM_ret = []
    for i in range(xs_sup.shape[0]):
        ysup = +m_abs_plot[i]
        yret = -m_abs_plot[i]
        lsup, = axM.plot(xs_sup[i], [ysup, ysup], color=color_supply, linewidth=1.0)
        lret, = axM.plot(xs_ret[i], [yret, yret], color=color_return, linewidth=1.0)
        linesM_sup.append(lsup)
        linesM_ret.append(lret)

    xh, yL, yD = _curve_at_t(t)
    base0 = np.zeros_like(xh)
    tot0 = yL + yD
    polyL = axH.fill_between(xh, base0, yL, color=color_loss, alpha=0.8)
    polyD = axH.fill_between(xh, yL, tot0, color=color_demand, alpha=0.8)
    lnH = axH.plot(xh, tot0, linewidth=1.0)[0]

    dt0 = pd.to_datetime(start_datetime) + pd.Timedelta(hours=int(t))
    Ts_min = float(np.nanmin(Ts)) if np.isfinite(np.nanmin(Ts)) else np.nan
    head_min = float(np.nanmin(_node_vec_attr(pressure_sup, t)) / (rho * g)) if np.isfinite(np.nanmin(_node_vec_attr(pressure_sup, t))) else np.nan
    Hinj = float(heat_inj[t]) if heat_inj is not None else np.nan
    Hloss = float(heat_loss[t]) if heat_loss is not None else np.nan
    head_max = float(np.nanmax(_node_vec_attr(pressure_sup, t)) / (rho * g)) if np.isfinite(np.nanmax(_node_vec_attr(pressure_sup, t))) else np.nan
    Ts_max = float(np.nanmax(Ts)) if np.isfinite(np.nanmax(Ts)) else np.nan
    if Hinj > 0:
        lossperc =  min (1, Hloss/Hinj)
    else:
        lossperc = 0
    Ppump = float(pump_pow[t]) if pump_pow is not None else np.nan

    _text(axInfo1, f"date/hour: {dt0}\nGround Temperature: {Tg_arr[t]:.2f}°C\n{name}")
    _text(axInfo2, f"lowest temperature at node: {Ts_min:.2f} °C\nlowest head (supply): {head_min:.2f} m")
    _text(axInfo3, f"pump power: {Ppump/1000:.2f} kW \nsupply station head = {head_max:.2f} m\nheat injected: {Hinj/1000:.2f} kW\n{100*lossperc:.2f}% due to losses\nsupply station temperature = {Ts_max:.2f} °C")

    def _update(_):
        t = int(sld.val)

        Ts = _node_vec_attr(temp_with_loss_sup, t)
        Tr = _node_vec_attr(temp_with_loss_ret, t)
        Ts_nl = _node_vec_attr(temp_no_loss_sup, t)
        Tr_nl = _node_vec_attr(temp_no_loss_ret, t)
        if Tg_arr.size:
            ground_line.set_ydata([Tg_arr[t], Tg_arr[t]])
            
        for i in range(xs_sup.shape[0]):
            linesT_sup[i].set_ydata([Ts[seg_s[i]], Ts[seg_e[i]]])
            linesT_ret[i].set_ydata([Tr[seg_e[i]], Tr[seg_s[i]]])
            linesTnl_sup[i].set_ydata([Ts_nl[seg_s[i]], Ts_nl[seg_e[i]]])
            linesTnl_ret[i].set_ydata([Tr_nl[seg_e[i]], Tr_nl[seg_s[i]]])

        scT_sup.set_offsets(np.c_[dist, Ts])
        scT_ret.set_offsets(np.c_[dist, Tr])
        scTnl_sup.set_offsets(np.c_[dist, Ts_nl])
        scTnl_ret.set_offsets(np.c_[dist, Tr_nl])

        Ps = _node_vec_attr(pressure_sup, t) / 100_000.0
        Pr = _node_vec_attr(pressure_ret, t) / 100_000.0

        for i in range(xs_sup.shape[0]):
            linesP_sup[i].set_ydata([Ps[seg_s[i]], Ps[seg_e[i]]])
            linesP_ret[i].set_ydata([Pr[seg_e[i]], Pr[seg_s[i]]])

        scP_sup.set_offsets(np.c_[dist, Ps])
        scP_ret.set_offsets(np.c_[dist, Pr])

        m_abs_plot = np.asarray([abs(_line_m_at_t(l, t)) for l in plot_lines], dtype=float)
        for i in range(xs_sup.shape[0]):
            ysup = +m_abs_plot[i]
            yret = -m_abs_plot[i]
            linesM_sup[i].set_ydata([ysup, ysup])
            linesM_ret[i].set_ydata([yret, yret])

        xh, yL, yD = _curve_at_t(t)
        tot = yL + yD

        nonlocal polyL, polyD
        polyL.remove()
        polyD.remove()
        polyL = axH.fill_between(xh, np.zeros_like(xh), yL, color=color_loss, alpha=0.8)
        polyD = axH.fill_between(xh, yL, tot, color=color_demand, alpha=0.8)
        lnH.set_data(xh, tot)

        dt0 = pd.to_datetime(start_datetime) + pd.Timedelta(hours=int(t))
        Ts_min = float(np.nanmin(Ts)) if np.isfinite(np.nanmin(Ts)) else np.nan
        Ts_max = float(np.nanmax(Ts)) if np.isfinite(np.nanmax(Ts)) else np.nan
        head_min = float(np.nanmin(_node_vec_attr(pressure_sup, t)) / (rho * g)) if np.isfinite(np.nanmin(_node_vec_attr(pressure_sup, t))) else np.nan
        head_max = float(np.nanmax(_node_vec_attr(pressure_sup, t)) / (rho * g)) if np.isfinite(np.nanmax(_node_vec_attr(pressure_sup, t))) else np.nan

        Hinj = float(heat_inj[t]) if heat_inj is not None else np.nan
        Hloss = float(heat_loss[t]) if heat_loss is not None else np.nan
        lossperc = min (1, Hloss/Hinj)
        Ppump = float(pump_pow[t]) if pump_pow is not None else np.nan

        _text(axInfo1, f"date/hour: {dt0}\nGround Temperature: {Tg_arr[t]:.2f}°C\n{name}")
        _text(axInfo2, f"lowest temperature at node: {Ts_min:.2f} °C\nlowest head (supply): {head_min:.2f} m")
        _text(axInfo3, f"pump power: {Ppump/1000:.2f} kW \nsupply station head = {head_max:.2f} m\nheat injected: {Hinj/1000:.2f} kW\n{100*lossperc:.2f}% due to losses\nsupply station temperature = {Ts_max:.2f} °C")

        fig.canvas.draw_idle()

    sld.on_changed(_update)
    plt.show()
    return fig
