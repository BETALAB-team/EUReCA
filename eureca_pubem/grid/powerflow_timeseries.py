from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional
import warnings
import numpy as np



from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np
import warnings

try:
    import numba as nb
    _HAS_NUMBA = True
except Exception:
    _HAS_NUMBA = False


@dataclass
class PowerFlowTSResults:
    V: np.ndarray
    P_inj: np.ndarray
    Q_inj: np.ndarray
    I_line: np.ndarray
    P_loss_line: np.ndarray
    E_loss_line_kWh: np.ndarray
    node_index: Dict[int, int]
    line_index: Dict[int, int]
    slack_bus: int
    converged: np.ndarray
    Vm_pu: Optional[np.ndarray] = None
    Va_rad: Optional[np.ndarray] = None


def q_from_p_pf(P, pf, lagging):
    P = np.asarray(P, dtype=float)
    if not (pf > 0 and pf <= 1):
        return np.zeros_like(P)
    if pf == 1:
        Q = np.zeros_like(P)
    else:
        Q = P * np.tan(np.arccos(pf))
    return Q if lagging else -Q


def compute_node_injections_timeseries(
    grid,
    pf_load_default=0.95,
    load_lagging_default=True,
    pf_pv_default=1.0,
    pv_lagging_default=False,
):
    T = None
    for n in getattr(grid, "nodes", []):
        if getattr(n, "node_type", "") == "building":
            bc = getattr(n, "building_consumption", None)
            bp = getattr(n, "building_production", None)
            if bc is not None and len(bc) > 0:
                T = len(bc); break
            if bp is not None and len(bp) > 0:
                T = len(bp); break
    if T is None:
        T = 1

    bus_ids = [int(getattr(n, "bus_id", 0)) for n in getattr(grid, "nodes", [])]
    node_index = {b: i for i, b in enumerate(bus_ids)}
    N = len(bus_ids)

    P_inj = np.zeros((T, N), dtype=float)
    Q_inj = np.zeros((T, N), dtype=float)

    for n in getattr(grid, "nodes", []):
        bid = int(getattr(n, "bus_id", 0))
        j = node_index.get(bid, None)
        if j is None:
            continue
        if getattr(n, "node_type", "") != "building":
            continue

        bc = getattr(n, "building_consumption", None)
        bp = getattr(n, "building_production", None)

        P_load = np.zeros(T) if bc is None else np.asarray(bc, dtype=float).reshape(-1)
        P_pv = np.zeros(T) if bp is None else np.asarray(bp, dtype=float).reshape(-1)

        if P_load.size != T:
            P_load = np.full(T, float(P_load[0])) if P_load.size == 1 else np.zeros(T)
        if P_pv.size != T:
            P_pv = np.full(T, float(P_pv[0])) if P_pv.size == 1 else np.zeros(T)

        P_load = np.maximum(P_load, 0.0)
        P_pv = np.maximum(P_pv, 0.0)

        pf_load = float(getattr(n, "pf_load", pf_load_default))
        load_lag = bool(getattr(n, "load_lagging", load_lagging_default))
        pf_pv = float(getattr(n, "pf_pv", pf_pv_default))
        pv_lag = bool(getattr(n, "pv_lagging", pv_lagging_default))

        Q_load = q_from_p_pf(P_load, pf_load, load_lag)

        inv = getattr(n, "inverter", None)
        if inv is None:
            P_gen = P_pv
        else:
            cap = float(getattr(inv, "capacity_w", 0.0))
            eff = float(getattr(inv, "efficiency", 1.0))
            if cap > 0 and np.isfinite(cap):
                P_gen = np.minimum(P_pv, cap) * eff
            else:
                P_gen = P_pv * eff

        Q_gen = np.zeros(T) if pf_pv == 1.0 else q_from_p_pf(P_gen, pf_pv, pv_lag)

        P_inj[:, j] = P_gen - P_load
        Q_inj[:, j] = Q_gen - Q_load

    return P_inj, Q_inj, node_index


def _pick_slack_bus(grid, slack_bus_id):
    supply = [int(getattr(n, "bus_id", 0)) for n in getattr(grid, "nodes", []) if getattr(n, "node_type", "") == "supply"]
    if slack_bus_id is not None and int(slack_bus_id) in supply:
        return int(slack_bus_id)
    if supply:
        return int(supply[0])
    nodes = getattr(grid, "nodes", [])
    return int(getattr(nodes[0], "bus_id", 0)) if nodes else 0


def _slack_setpoint(grid, slack_bus):
    V = 1.0
    ang = 0.0
    for n in getattr(grid, "nodes", []):
        if int(getattr(n, "bus_id", 0)) == int(slack_bus):
            sv = getattr(n, "set_voltage", None)
            sa = getattr(n, "set_voltage_angle", None)
            if sv is not None and np.isfinite(float(sv)) and float(sv) > 0:
                V = float(sv)
            if sa is not None and np.isfinite(float(sa)):
                ang = float(sa)
            break
    return V, ang


# -------------------------
# Array-compiled grid
# -------------------------
@dataclass
class CompiledRadialGrid:
    # node mapping
    bus_ids: np.ndarray             # (N,)
    node_index: Dict[int, int]      # bus_id -> idx
    slack_bus: int
    slack_idx: int

    # tree structure
    parent_idx: np.ndarray          # (N,) int32, -1 for slack/unreached
    parent_line_idx: np.ndarray     # (N,) int32, -1 for slack/unreached
    order: np.ndarray              # (N_reached,) int32 node indices
    post: np.ndarray               # (N_reached,) int32 node indices
    child_ptr: np.ndarray          # (N+1,) int32
    child_idx: np.ndarray          # (N_reached-1,) int32 (children of each node)

    # lines (all L)
    from_idx: np.ndarray           # (L,) int32
    to_idx: np.ndarray             # (L,) int32
    z_line: np.ndarray             # (L,) complex128
    r_line: np.ndarray             # (L,) float64
    line_id: np.ndarray            # (L,) int32
    line_index: Dict[int, int]     # line_id -> li

    # for forward sweep we need impedance of each child-edge by node:
    z_to_parent: np.ndarray        # (N,) complex128; z of parent edge (0 for slack/unreached)


def compile_radial_grid(grid, slack_bus_id=None) -> CompiledRadialGrid:
    nodes = list(getattr(grid, "nodes", []))
    lines = list(getattr(grid, "lines", []))

    bus_ids = np.asarray([int(getattr(n, "bus_id", 0)) for n in nodes], dtype=np.int64)
    node_index = {int(b): i for i, b in enumerate(bus_ids)}
    N = int(bus_ids.size)

    slack_bus = _pick_slack_bus(grid, slack_bus_id)
    slack_idx = int(node_index.get(int(slack_bus), 0))

    # line arrays
    L = len(lines)
    from_bus = np.empty(L, dtype=np.int64)
    to_bus = np.empty(L, dtype=np.int64)
    line_id = np.empty(L, dtype=np.int64)
    z_line = np.empty(L, dtype=np.complex128)
    r_line = np.empty(L, dtype=np.float64)

    for i, l in enumerate(lines):
        u = int(getattr(l, "from_bus", 0))
        v = int(getattr(l, "to_bus", 0))
        lid = int(getattr(l, "line_id", 0))
        rpkm = float(getattr(l, "r_ohm_per_km", 0.0) or 0.0)
        xpkm = float(getattr(l, "x_ohm_per_km", 0.0) or 0.0)
        length_km = float(getattr(l, "length_m", 0.0) or 0.0) / 1000.0
        z = (rpkm + 1j * xpkm) * length_km
        if not (np.isfinite(z.real) and np.isfinite(z.imag)) or abs(z) <= 0:
            z = 1e9 + 0j

        from_bus[i] = u
        to_bus[i] = v
        line_id[i] = lid
        z_line[i] = z
        r_line[i] = float(np.real(z))

    line_index = {int(lid): i for i, lid in enumerate(line_id.tolist())}

    # adjacency using indices for BFS tree build
    # build edge lists (undirected) in index space
    uu = []
    vv = []
    ll = []
    for i in range(L):
        u = int(from_bus[i]); v = int(to_bus[i])
        if u not in node_index or v not in node_index:
            continue
        iu = int(node_index[u]); iv = int(node_index[v])
        uu.append(iu); vv.append(iv); ll.append(i)

    uu = np.asarray(uu, dtype=np.int32)
    vv = np.asarray(vv, dtype=np.int32)
    ll = np.asarray(ll, dtype=np.int32)
    M = int(uu.size)

    # build CSR adjacency for BFS
    deg = np.zeros(N, dtype=np.int32)
    for e in range(M):
        deg[uu[e]] += 1
        deg[vv[e]] += 1
    ptr = np.zeros(N + 1, dtype=np.int32)
    np.cumsum(deg, out=ptr[1:])
    adj_nbr = np.empty(ptr[-1], dtype=np.int32)
    adj_edge = np.empty(ptr[-1], dtype=np.int32)
    cur = ptr[:-1].copy()
    for e in range(M):
        a = uu[e]; b = vv[e]; li = ll[e]
        pa = cur[a]; adj_nbr[pa] = b; adj_edge[pa] = li; cur[a] += 1
        pb = cur[b]; adj_nbr[pb] = a; adj_edge[pb] = li; cur[b] += 1

    # BFS tree
    parent_idx = np.full(N, -1, dtype=np.int32)
    parent_line_idx = np.full(N, -1, dtype=np.int32)
    order = np.empty(N, dtype=np.int32)
    q = np.empty(N, dtype=np.int32)

    parent_idx[slack_idx] = slack_idx
    parent_line_idx[slack_idx] = -1
    qh = 0; qt = 0
    q[qt] = slack_idx; qt += 1
    oc = 0

    while qh < qt:
        u = q[qh]; qh += 1
        order[oc] = u; oc += 1
        for p in range(ptr[u], ptr[u + 1]):
            v = adj_nbr[p]
            if parent_idx[v] != -1:
                continue
            parent_idx[v] = u
            parent_line_idx[v] = adj_edge[p]
            q[qt] = v; qt += 1

    order = order[:oc].copy()
    post = order[::-1].copy()

    # children CSR
    # count children only for reached nodes
    child_count = np.zeros(N, dtype=np.int32)
    for v in order:
        if v == slack_idx:
            continue
        p = parent_idx[v]
        if p >= 0 and p != v:
            child_count[p] += 1

    child_ptr = np.zeros(N + 1, dtype=np.int32)
    np.cumsum(child_count, out=child_ptr[1:])
    child_idx = np.empty(max(int(child_ptr[-1]), 0), dtype=np.int32)
    curc = child_ptr[:-1].copy()
    for v in order:
        if v == slack_idx:
            continue
        p = parent_idx[v]
        if p >= 0 and p != v:
            pos = curc[p]
            child_idx[pos] = v
            curc[p] += 1

    # z_to_parent by node
    z_to_parent = np.zeros(N, dtype=np.complex128)
    for v in order:
        if v == slack_idx:
            continue
        li = int(parent_line_idx[v])
        if li >= 0:
            z_to_parent[v] = z_line[li]
        else:
            z_to_parent[v] = 1e9 + 0j

    # line endpoint indices (all lines; directed as stored)
    from_idx = np.full(L, -1, dtype=np.int32)
    to_idx = np.full(L, -1, dtype=np.int32)
    for i in range(L):
        u = int(from_bus[i]); v = int(to_bus[i])
        from_idx[i] = int(node_index.get(u, -1))
        to_idx[i] = int(node_index.get(v, -1))

    return CompiledRadialGrid(
        bus_ids=bus_ids.astype(np.int64),
        node_index=node_index,
        slack_bus=int(slack_bus),
        slack_idx=int(slack_idx),
        parent_idx=parent_idx,
        parent_line_idx=parent_line_idx,
        order=order,
        post=post,
        child_ptr=child_ptr,
        child_idx=child_idx,
        from_idx=from_idx,
        to_idx=to_idx,
        z_line=z_line,
        r_line=r_line,
        line_id=line_id.astype(np.int32),
        line_index=line_index,
        z_to_parent=z_to_parent,
    )


if _HAS_NUMBA:
    @nb.njit(cache=False, fastmath=False)
    def _bfs_solve_timeseries_numba(
        P_inj, Q_inj,
        slack_idx, Vslack_V, ang_slack,
        order, post, child_ptr, child_idx, parent_idx,
        z_to_parent,
        from_idx, to_idx, z_line, r_line,
        timestep_hours, max_iter, tol
    ):
        T, N = P_inj.shape
        L = from_idx.size

        Vmag = np.full((T, N), np.nan, dtype=np.float64)
        Vm_pu = np.full((T, N), np.nan, dtype=np.float64)
        Va = np.full((T, N), np.nan, dtype=np.float64)
        I_line = np.full((T, L), np.nan, dtype=np.float64)
        P_loss = np.full((T, L), np.nan, dtype=np.float64)
        E_loss = np.full((T, L), np.nan, dtype=np.float64)
        converged = np.zeros(T, dtype=np.bool_)

        V = np.empty(N, dtype=np.complex128)
        Vprev = np.empty(N, dtype=np.complex128)
        I_inj = np.zeros(N, dtype=np.complex128)
        I_branch = np.zeros(N, dtype=np.complex128)

        Vslack = Vslack_V * (np.cos(ang_slack) + 1j*np.sin(ang_slack))

        for t in range(T):
            # init V
            for i in range(N):
                V[i] = Vslack
            V[slack_idx] = Vslack

            ok = False
            for it in range(max_iter):
                for i in range(N):
                    Vprev[i] = V[i]

                # injections current: I = -conj(S/V)
                for i in range(N):
                    I_inj[i] = 0.0 + 0.0j
                for i in range(N):
                    if i == slack_idx:
                        continue
                    Vr = V[i].real
                    Vi = V[i].imag
                    denom = Vr*Vr + Vi*Vi
                    if denom < 1e-12:
                        denom = 1e-12
                    # S = P + jQ
                    P = P_inj[t, i]
                    Q = Q_inj[t, i]
                    # S/V = (P+jQ)/(Vr+jVi) = ( (P+jQ)*(Vr-jVi) )/denom
                    a = (P*Vr + Q*Vi) / denom
                    b = (Q*Vr - P*Vi) / denom
                    # conj(S/V) = a - jb
                    I_inj[i] = -(a - 1j*b)

                # backward sweep (post-order)
                for idx in range(post.size):
                    b = post[idx]
                    Itot = I_inj[b]
                    c0 = child_ptr[b]
                    c1 = child_ptr[b+1]
                    for p in range(c0, c1):
                        c = child_idx[p]
                        Itot += I_branch[c]
                    I_branch[b] = Itot

                # forward sweep (order)
                for idx in range(order.size):
                    b = order[idx]
                    if b == slack_idx:
                        continue
                    p = parent_idx[b]
                    if p < 0:
                        continue
                    V[b] = V[p] - z_to_parent[b]*I_branch[b]

                # convergence
                maxdv = 0.0
                for i in range(N):
                    dv = V[i] - Vprev[i]
                    ad = (dv.real*dv.real + dv.imag*dv.imag) ** 0.5
                    if ad > maxdv:
                        maxdv = ad
                if maxdv < tol:
                    ok = True
                    break

            converged[t] = ok
            # store node voltages
            for i in range(N):
                vm = (V[i].real*V[i].real + V[i].imag*V[i].imag) ** 0.5
                Vmag[t, i] = vm
                Vm_pu[t, i] = vm / (Vslack_V if Vslack_V > 1e-12 else 1e-12)
                Va[t, i] = np.arctan2(V[i].imag, V[i].real)

            # line currents + losses (using line direction as stored)
            for li in range(L):
                u = from_idx[li]
                v = to_idx[li]
                if u < 0 or v < 0:
                    continue
                z = z_line[li]
                # I = (Vu - Vv)/z
                dV = V[u] - V[v]
                denom = z.real*z.real + z.imag*z.imag
                if denom < 1e-24:
                    continue
                Iuv = dV / z
                Iabs = (Iuv.real*Iuv.real + Iuv.imag*Iuv.imag) ** 0.5
                I_line[t, li] = Iabs
                Pl = Iabs*Iabs * r_line[li]
                P_loss[t, li] = Pl
                E_loss[t, li] = (Pl * timestep_hours) / 1000.0

        return Vmag, Vm_pu, Va, I_line, P_loss, E_loss, converged


def solve_grid_timeseries_bfs_fast(
    grid,
    slack_bus_id=None,
    timestep_hours=1.0,
    pf_load_default=0.95,
    load_lagging_default=True,
    pf_pv_default=1.0,
    pv_lagging_default=False,
    max_iter=1000,
    tol=1e-6,
):
    if not _HAS_NUMBA:
        raise ImportError("numba is not installed. Install it with: pip install numba")

    slack_bus = _pick_slack_bus(grid, slack_bus_id)
    Vslack_V, ang_slack = _slack_setpoint(grid, slack_bus)

    P_inj, Q_inj, node_index = compute_node_injections_timeseries(
        grid,
        pf_load_default=pf_load_default,
        load_lagging_default=load_lagging_default,
        pf_pv_default=pf_pv_default,
        pv_lagging_default=pv_lagging_default,
    )

    cg = compile_radial_grid(grid, slack_bus_id=slack_bus)
    Vmag, Vm_pu, Va, I_line, P_loss_line, E_loss_line_kWh, converged = _bfs_solve_timeseries_numba(
        P_inj.astype(np.float64),
        Q_inj.astype(np.float64),
        cg.slack_idx,
        float(Vslack_V),
        float(ang_slack),
        cg.order,
        cg.post,
        cg.child_ptr,
        cg.child_idx,
        cg.parent_idx,
        cg.z_to_parent,
        cg.from_idx,
        cg.to_idx,
        cg.z_line,
        cg.r_line,
        float(timestep_hours),
        int(max_iter),
        float(tol),
    )

    res = PowerFlowTSResults(
        V=Vmag,
        P_inj=P_inj,
        Q_inj=Q_inj,
        I_line=I_line,
        P_loss_line=P_loss_line,
        E_loss_line_kWh=E_loss_line_kWh,
        node_index=node_index,
        line_index=cg.line_index,
        slack_bus=int(slack_bus),
        converged=converged,
        Vm_pu=Vm_pu,
        Va_rad=Va,
    )

    # supply metrics
    res.P_supply = -P_inj.sum(axis=1, keepdims=True) + P_loss_line.sum(axis=1, keepdims=True)
    res.Q_supply = -Q_inj.sum(axis=1, keepdims=True)
    res.S_supply = np.sqrt(res.P_supply**2 + res.Q_supply**2)

    if not np.all(converged):
        # warn once, not 8760 times
        bad = np.where(~converged)[0]

    return res


@dataclass
class PowerFlowTSResults:
    V: np.ndarray
    P_inj: np.ndarray
    Q_inj: np.ndarray
    I_line: np.ndarray
    P_loss_line: np.ndarray
    E_loss_line_kWh: np.ndarray
    node_index: Dict[int, int]
    line_index: Dict[int, int]
    slack_bus: int
    converged: np.ndarray
    Vm_pu: Optional[np.ndarray] = None
    Va_rad: Optional[np.ndarray] = None


def q_from_p_pf(P, pf, lagging):
    P = np.asarray(P, dtype=float)
    if not (pf > 0 and pf <= 1):
        return np.zeros_like(P)
    if pf == 1:
        Q = np.zeros_like(P)
    else:
        Q = P * np.tan(np.arccos(pf))
    return Q if lagging else -Q


def compute_node_injections_timeseries(
    grid,
    pf_load_default=0.95,
    load_lagging_default=True,
    pf_pv_default=1.0,
    pv_lagging_default=False,
):
    T = None
    for n in getattr(grid, "nodes", []):
        if getattr(n, "node_type", "") == "building":
            bc = getattr(n, "building_consumption", None)
            bp = getattr(n, "building_production", None)
            if bc is not None and len(bc) > 0:
                T = len(bc)
                break
            if bp is not None and len(bp) > 0:
                T = len(bp)
                break
    if T is None:
        T = 1

    bus_ids = [int(getattr(n, "bus_id", 0)) for n in getattr(grid, "nodes", [])]
    node_index = {b: i for i, b in enumerate(bus_ids)}
    N = len(bus_ids)

    P_inj = np.zeros((T, N), dtype=float)
    Q_inj = np.zeros((T, N), dtype=float)

    for n in getattr(grid, "nodes", []):
        bid = int(getattr(n, "bus_id", 0))
        j = node_index.get(bid, None)
        if j is None:
            continue
        if getattr(n, "node_type", "") != "building":
            continue

        bc = getattr(n, "building_consumption", None)
        bp = getattr(n, "building_production", None)

        P_load = np.zeros(T) if bc is None else np.asarray(bc, dtype=float).reshape(-1)
        P_pv = np.zeros(T) if bp is None else np.asarray(bp, dtype=float).reshape(-1)

        if P_load.size != T:
            P_load = np.full(T, float(P_load[0])) if P_load.size == 1 else np.zeros(T)
        if P_pv.size != T:
            P_pv = np.full(T, float(P_pv[0])) if P_pv.size == 1 else np.zeros(T)

        P_load = np.maximum(P_load, 0.0)
        P_pv = np.maximum(P_pv, 0.0)

        pf_load = float(getattr(n, "pf_load", pf_load_default))
        load_lag = bool(getattr(n, "load_lagging", load_lagging_default))
        pf_pv = float(getattr(n, "pf_pv", pf_pv_default))
        pv_lag = bool(getattr(n, "pv_lagging", pv_lagging_default))

        Q_load = q_from_p_pf(P_load, pf_load, load_lag)

        inv = getattr(n, "inverter", None)
        if inv is None:
            P_gen = P_pv
        else:
            cap = float(getattr(inv, "capacity_w", 0.0))
            eff = float(getattr(inv, "efficiency", 1.0))
            if cap > 0 and np.isfinite(cap):
                P_gen = np.minimum(P_pv, cap) * eff
            else:
                P_gen = P_pv * eff

        if pf_pv == 1.0:
            Q_gen = np.zeros(T)
        else:
            Q_gen = q_from_p_pf(P_gen, pf_pv, pv_lag)

        P_inj[:, j] = P_gen - P_load
        Q_inj[:, j] = Q_gen - Q_load

    return P_inj, Q_inj, node_index


def _pick_slack_bus(grid, slack_bus_id):
    supply = [int(getattr(n, "bus_id", 0)) for n in getattr(grid, "nodes", []) if getattr(n, "node_type", "") == "supply"]
    if slack_bus_id is not None and int(slack_bus_id) in supply:
        return int(slack_bus_id)
    if supply:
        return int(supply[0])
    nodes = getattr(grid, "nodes", [])
    return int(getattr(nodes[0], "bus_id", 0)) if nodes else 0


def _slack_setpoint(grid, slack_bus):
    V = 400.0
    ang = 0.0
    for n in getattr(grid, "nodes", []):
        if int(getattr(n, "bus_id", 0)) == int(slack_bus):
            sv = getattr(n, "set_voltage", None)
            sa = getattr(n, "set_voltage_angle", None)
            if sv is not None and np.isfinite(float(sv)) and float(sv) > 0:
                V = float(sv)
            if sa is not None and np.isfinite(float(sa)):
                ang = float(sa)
            break
    return V, ang


def _build_adjacency(lines):
    adj = {}
    for l in lines:
        u = int(getattr(l, "from_bus", 0))
        v = int(getattr(l, "to_bus", 0))
        lid = int(getattr(l, "line_id", 0))
        adj.setdefault(u, []).append((v, lid))
        adj.setdefault(v, []).append((u, lid))
    return adj


def make_grid_with_lines_removed(grid, drop_line_ids: List[int]):
    drop = set(int(x) for x in (drop_line_ids or []))
    kept = [l for l in getattr(grid, "lines", []) if int(getattr(l, "line_id", 0)) not in drop]
    return type(grid)(nodes=getattr(grid, "nodes", []), lines=kept)


def _build_ybus_ohm(grid, node_index):
    N = len(node_index)
    Y = np.zeros((N, N), dtype=complex)
    for l in getattr(grid, "lines", []):
        i = node_index.get(int(getattr(l, "from_bus", 0)), None)
        k = node_index.get(int(getattr(l, "to_bus", 0)), None)
        if i is None or k is None:
            continue
        rpkm = float(getattr(l, "r_ohm_per_km", 0.0) or 0.0)
        xpkm = float(getattr(l, "x_ohm_per_km", 0.0) or 0.0)
        z = (rpkm + 1j * xpkm) * (float(getattr(l, "length_m", 0.0)) / 1000.0)
        if not (np.isfinite(z.real) and np.isfinite(z.imag)) or abs(z) <= 0:
            continue
        y = 1.0 / z
        Y[i, i] += y
        Y[k, k] += y
        Y[i, k] -= y
        Y[k, i] -= y
    return Y


def _bfs_tree_edges(grid, slack_bus):
    lines = getattr(grid, "lines", [])
    adj = _build_adjacency(lines)
    buses = [int(getattr(n, "bus_id", 0)) for n in getattr(grid, "nodes", [])]
    bus_set = set(buses)

    parent_of = {slack_bus: None}
    parent_line_of = {slack_bus: None}
    order = [slack_bus]
    q = [slack_bus]
    seen = {slack_bus}

    while q:
        u = q.pop(0)
        for v, lid in adj.get(u, []):
            if v not in bus_set:
                continue
            if v in seen:
                continue
            seen.add(v)
            parent_of[v] = u
            parent_line_of[v] = lid
            q.append(v)
            order.append(v)

    children = {b: [] for b in bus_set}
    for b, p in parent_of.items():
        if p is not None:
            children[p].append(b)

    post = list(reversed(order))
    return parent_of, parent_line_of, order, post, children


def solve_grid_timeseries_bfs(
    grid,
    slack_bus_id=None,
    timestep_hours=1.0,
    pf_load_default=0.95,
    load_lagging_default=True,
    pf_pv_default=1.0,
    pv_lagging_default=False,
    max_iter=30,
    tol=1e-6,
):
    slack_bus = _pick_slack_bus(grid, slack_bus_id)
    Vslack_V, ang_slack = _slack_setpoint(grid, slack_bus)

    P_inj, Q_inj, node_index = compute_node_injections_timeseries(
        grid,
        pf_load_default=pf_load_default,
        load_lagging_default=load_lagging_default,
        pf_pv_default=pf_pv_default,
        pv_lagging_default=pv_lagging_default,
    )

    T, N = P_inj.shape
    lines = list(getattr(grid, "lines", []))
    L = len(lines)
    line_index = {int(getattr(l, "line_id", 0)): i for i, l in enumerate(lines)}

    parent_of, parent_line_of, order, post, children = _bfs_tree_edges(grid, slack_bus)

    z_by_lid = {}
    for l in lines:
        lid = int(getattr(l, "line_id", 0))
        rpkm = float(getattr(l, "r_ohm_per_km", 0.0) or 0.0)
        xpkm = float(getattr(l, "x_ohm_per_km", 0.0) or 0.0)
        z = (rpkm + 1j * xpkm) * (float(getattr(l, "length_m", 0.0)) / 1000.0)
        if abs(z) <= 0:
            z = 1e9 + 0j
        z_by_lid[lid] = z

    converged = np.zeros(T, dtype=bool)
    Vmag = np.full((T, N), np.nan, dtype=float)
    I_line = np.full((T, L), np.nan, dtype=float)
    P_loss_line = np.full((T, L), np.nan, dtype=float)
    E_loss_line_kWh = np.full((T, L), np.nan, dtype=float)
    Vm_pu = np.full((T, N), np.nan, dtype=float)
    Va_rad = np.full((T, N), np.nan, dtype=float)

    for t in range(T):
        S = (P_inj[t, :] + 1j * Q_inj[t, :]).astype(complex)
        V = np.ones(N, dtype=complex) * (Vslack_V * np.exp(1j * ang_slack))
        V[node_index.get(slack_bus, 0)] = Vslack_V * np.exp(1j * ang_slack)

        ok = False
        for _ in range(max_iter):
            V_prev = V.copy()

            I_inj = np.zeros(N, dtype=complex)
            mask = np.ones(N, dtype=bool)
            mask[node_index.get(slack_bus, 0)] = False
            V_use = V.copy()
            V_use[np.abs(V_use) < 1e-6] = 1e-6 + 0j
            I_inj[mask] = -np.conj(S[mask] / V_use[mask])

            I_branch = {b: 0j for b in parent_of.keys()}
            for b in post:
                idxb = node_index.get(b, None)
                if idxb is None:
                    continue
                I_total = I_inj[idxb]
                for c in children.get(b, []):
                    I_total += I_branch.get(c, 0j)
                I_branch[b] = I_total

            for b in order:
                if b == slack_bus:
                    continue
                p = parent_of.get(b, None)
                lid = parent_line_of.get(b, None)
                if p is None or lid is None:
                    continue
                ip = node_index.get(p, None)
                ib = node_index.get(b, None)
                if ip is None or ib is None:
                    continue
                z = z_by_lid.get(int(lid), 1e9 + 0j)
                V[ib] = V[ip] - z * I_branch.get(b, 0j)

            if np.nanmax(np.abs(V - V_prev)) < tol:
                ok = True
                break

        converged[t] = ok
        Vmag[t, :] = np.abs(V)
        Vm_pu[t, :] = np.abs(V) / max(Vslack_V, 1e-9)
        Va_rad[t, :] = np.angle(V)

        for li, l in enumerate(lines):
            u = node_index.get(int(getattr(l, "from_bus", 0)), None)
            v = node_index.get(int(getattr(l, "to_bus", 0)), None)
            if u is None or v is None:
                continue
            z = z_by_lid.get(int(getattr(l, "line_id", 0)), 1e9 + 0j)
            if abs(z) <= 0:
                continue
            Iuv = (V[u] - V[v]) / z
            Iabs = float(np.abs(Iuv))
            I_line[t, li] = Iabs
            R = float(np.real(z))
            Ploss = (Iabs * Iabs) * R
            P_loss_line[t, li] = Ploss
            E_loss_line_kWh[t, li] = (Ploss * timestep_hours) / 1000.0

        if not ok:
            warnings.warn(f"BFS did not converge at t={t}")

    res = PowerFlowTSResults(
        V=Vmag,
        P_inj=P_inj,
        Q_inj=Q_inj,
        I_line=I_line,
        P_loss_line=P_loss_line,
        E_loss_line_kWh=E_loss_line_kWh,
        node_index=node_index,
        line_index=line_index,
        slack_bus=slack_bus,
        converged=converged,
        Vm_pu=Vm_pu,
        Va_rad=Va_rad,
    )
    res.P_supply = -P_inj.sum(axis=1, keepdims=True)+P_loss_line.sum(axis=1, keepdims=True) 
    res.Q_supply = -Q_inj.sum(axis=1, keepdims=True)
    res.S_supply = np.sqrt(res.P_supply**2 + res.Q_supply**2)
    
    return res

def apply_results_to_grid(grid, res, *, prefix=""):
    ni = getattr(res, "node_index", {})
    li = getattr(res, "line_index", {})

    V = getattr(res, "V", None)
    P_inj = getattr(res, "P_inj", None)
    Q_inj = getattr(res, "Q_inj", None)

    I_line = getattr(res, "I_line", None)
    P_loss_line = getattr(res, "P_loss_line", None)
    E_loss_line_kWh = getattr(res, "E_loss_line_kWh", None)

    for n in getattr(grid, "nodes", []):
        bid = int(getattr(n, "bus_id", 0))
        j = ni.get(bid, None)
        if j is None:
            continue

        if V is not None:
            object.__setattr__(n, f"{prefix}V", V[:, j])
        if P_inj is not None:
            object.__setattr__(n, f"{prefix}P_inj", P_inj[:, j])
        if Q_inj is not None:
            object.__setattr__(n, f"{prefix}Q_inj", Q_inj[:, j])

    for l in getattr(grid, "lines", []):
        lid = int(getattr(l, "line_id", 0))
        k = li.get(lid, None)
        if k is None:
            continue

        if I_line is not None:
            setattr(l, f"{prefix}I", I_line[:, k])
        if P_loss_line is not None:
            setattr(l, f"{prefix}P_loss", P_loss_line[:, k])
        if E_loss_line_kWh is not None:
            setattr(l, f"{prefix}E_loss_kWh", E_loss_line_kWh[:, k])

    setattr(grid, f"{prefix}results", res)
    return grid

def _build_ybus_pu(grid, node_index, V_base, S_base):
    Y_ohm = _build_ybus_ohm(grid, node_index)
    Z_base = (V_base * V_base) / max(S_base, 1e-12)
    Y_pu = Y_ohm * Z_base
    return Y_pu


def _nr_solve_pq_fast(
    Y_pu,
    slack_idx,
    Vslack_pu,
    P_spec_pu,
    Q_spec_pu,
    pq,
    G,
    B,
    max_iter,
    tol,
    Vm0=None,
    Va0=None,
    rebuild_J_every=1,
):
    N = Y_pu.shape[0]
    npq = pq.size

    if Vm0 is None or Va0 is None:
        Vm = np.ones(N, dtype=float)
        Va = np.zeros(N, dtype=float)
    else:
        Vm = np.asarray(Vm0, dtype=float).copy()
        Va = np.asarray(Va0, dtype=float).copy()
        if Vm.size != N or Va.size != N:
            Vm = np.ones(N, dtype=float)
            Va = np.zeros(N, dtype=float)

    Vm[slack_idx] = np.abs(Vslack_pu)
    Va[slack_idx] = np.angle(Vslack_pu)

    J = None

    for it in range(max_iter):
        V = Vm * np.exp(1j * Va)
        I = Y_pu @ V
        S = V * np.conj(I)
        P = S.real
        Q = S.imag

        dP = P_spec_pu - P
        dQ = Q_spec_pu - Q
        mis = np.concatenate([dP[pq], dQ[pq]])
        mmax = float(np.max(np.abs(mis))) if mis.size else 0.0
        if mmax < tol:
            return Vm, Va, True, it + 1

        if (J is None) or (rebuild_J_every <= 1) or (it % rebuild_J_every == 0):
            npq = pq.size
            H = np.zeros((npq, npq), dtype=float)
            Nmat = np.zeros((npq, npq), dtype=float)
            M = np.zeros((npq, npq), dtype=float)
            L = np.zeros((npq, npq), dtype=float)

            for a in range(npq):
                i = int(pq[a])
                for b in range(npq):
                    k = int(pq[b])
                    if i == k:
                        H[a, b] = -Q[i] - (B[i, i] * Vm[i] * Vm[i])
                        Nmat[a, b] = P[i] / max(Vm[i], 1e-12) + (G[i, i] * Vm[i])
                        M[a, b] = P[i] - (G[i, i] * Vm[i] * Vm[i])
                        L[a, b] = Q[i] / max(Vm[i], 1e-12) - (B[i, i] * Vm[i])
                    else:
                        th = Va[i] - Va[k]
                        c = np.cos(th)
                        s = np.sin(th)
                        H[a, b] = Vm[i] * Vm[k] * (G[i, k] * s - B[i, k] * c)
                        Nmat[a, b] = Vm[i] * (G[i, k] * c + B[i, k] * s)
                        M[a, b] = -Vm[i] * Vm[k] * (G[i, k] * c + B[i, k] * s)
                        L[a, b] = Vm[i] * (G[i, k] * s - B[i, k] * c)

            J = np.block([[H, Nmat], [M, L]])

        try:
            dx = np.linalg.solve(J, mis)
        except Exception:
            return Vm, Va, False, it + 1

        dth = dx[:pq.size]
        dV = dx[pq.size:]

        Va[pq] = Va[pq] + dth
        Vm[pq] = Vm[pq] + dV

        bad = (~np.isfinite(Vm)) | (Vm <= 1e-6)
        if np.any(bad):
            Vm[bad] = 1e-6

        Vm[slack_idx] = np.abs(Vslack_pu)
        Va[slack_idx] = np.angle(Vslack_pu)

    return Vm, Va, False, max_iter


def solve_grid_timeseries_newton_raphson(
    grid,
    slack_bus_id=None,
    timestep_hours=1.0,
    pf_load_default=0.95,
    load_lagging_default=True,
    pf_pv_default=1.0,
    pv_lagging_default=False,
    S_base_va=1_000_000.0,
    tol_pu=1e-3,
    max_iter=25,
    rebuild_J_every=2,
    Vm0=None,
    Va0=None,
):
    slack_bus = _pick_slack_bus(grid, slack_bus_id)
    V_slack_V, ang_slack = _slack_setpoint(grid, slack_bus)

    P_inj, Q_inj, node_index = compute_node_injections_timeseries(
        grid,
        pf_load_default=pf_load_default,
        load_lagging_default=load_lagging_default,
        pf_pv_default=pf_pv_default,
        pv_lagging_default=pv_lagging_default,
    )

    T, N = P_inj.shape
    lines = list(getattr(grid, "lines", []))
    L = len(lines)
    line_index = {int(getattr(l, "line_id", 0)): i for i, l in enumerate(lines)}

    V_base = float(V_slack_V) if np.isfinite(V_slack_V) and V_slack_V > 0 else 1.0
    S_base = float(S_base_va) if np.isfinite(S_base_va) and S_base_va > 0 else 1_000_000.0

    Y_pu = _build_ybus_pu(grid, node_index, V_base, S_base)
    slack_idx = node_index.get(int(slack_bus), 0)
    Vslack_pu = (V_slack_V / max(V_base, 1e-12)) * np.exp(1j * ang_slack)

    G = Y_pu.real
    B = Y_pu.imag
    pq = np.array([i for i in range(N) if i != slack_idx], dtype=int)

    Vmag_V = np.full((T, N), np.nan, dtype=float)
    I_line = np.full((T, L), np.nan, dtype=float)
    P_loss_line = np.full((T, L), np.nan, dtype=float)
    E_loss_line_kWh = np.full((T, L), np.nan, dtype=float)
    converged = np.zeros(T, dtype=bool)
    Vm_hist = np.full((T, N), np.nan, dtype=float)
    Va_hist = np.full((T, N), np.nan, dtype=float)

    Vm_guess = None if Vm0 is None else np.asarray(Vm0, dtype=float).copy()
    Va_guess = None if Va0 is None else np.asarray(Va0, dtype=float).copy()

    for t in range(T):
        P_spec_pu = (P_inj[t, :] / S_base).astype(float)
        Q_spec_pu = (Q_inj[t, :] / S_base).astype(float)

        Vm, Va, ok, _ = _nr_solve_pq_fast(
            Y_pu=Y_pu,
            slack_idx=slack_idx,
            Vslack_pu=Vslack_pu,
            P_spec_pu=P_spec_pu,
            Q_spec_pu=Q_spec_pu,
            pq=pq,
            G=G,
            B=B,
            max_iter=max_iter,
            tol=tol_pu,
            Vm0=Vm_guess,
            Va0=Va_guess,
            rebuild_J_every=rebuild_J_every,
        )

        converged[t] = bool(ok)
        if ok:
            Vm_guess = Vm
            Va_guess = Va
        else:
            Vm_guess = None
            Va_guess = None
            # warnings.warn(f"NR did not converge at t={t}")

        Vpu = Vm * np.exp(1j * Va)
        Vcplx = Vpu * V_base
        Vmag_V[t, :] = np.abs(Vcplx)
        Vm_hist[t, :] = Vm
        Va_hist[t, :] = Va

        for li, l in enumerate(lines):
            i = node_index.get(int(getattr(l, "from_bus", 0)), None)
            k = node_index.get(int(getattr(l, "to_bus", 0)), None)
            if i is None or k is None:
                continue

            rpkm = float(getattr(l, "r_ohm_per_km", 0.0) or 0.0)
            xpkm = float(getattr(l, "x_ohm_per_km", 0.0) or 0.0)
            z = (rpkm + 1j * xpkm) * (float(getattr(l, "length_m", 0.0)) / 1000.0)
            if abs(z) <= 0:
                continue
            y = 1.0 / z

            Iik = (Vcplx[i] - Vcplx[k]) * y
            Iabs = float(np.abs(Iik))
            I_line[t, li] = Iabs

            R = float(np.real(z))
            Ploss = (Iabs * Iabs) * R
            P_loss_line[t, li] = Ploss
            E_loss_line_kWh[t, li] = (Ploss * timestep_hours) / 1000.0

    return PowerFlowTSResults(
        V=Vmag_V,
        P_inj=P_inj,
        Q_inj=Q_inj,
        I_line=I_line,
        P_loss_line=P_loss_line,
        E_loss_line_kWh=E_loss_line_kWh,
        node_index=node_index,
        line_index=line_index,
        slack_bus=slack_bus,
        converged=converged,
        Vm_pu=Vm_hist,
        Va_rad=Va_hist,
    )


def solve_grid_timeseries_auto(
    grid,
    slack_bus_id=None,
    timestep_hours=1.0,
    pf_load_default=0.95,
    load_lagging_default=True,
    pf_pv_default=1.0,
    pv_lagging_default=False,
    bfs_max_iter=30,
    bfs_tol=1e-6,
    nr_S_base_va=1_000_000.0,
    nr_tol_pu=1e-8,
    nr_max_iter=20,
    nr_rebuild_J_every=2,
):
    kind_obj = getattr(grid, "topology", None)
    kind = getattr(kind_obj, "kind", None)
    drop_line_ids = getattr(kind_obj, "drop_line_ids", None)
    slack = getattr(kind_obj, "slack_bus", None)

    if slack_bus_id is None and slack is not None:
        slack_bus_id = int(slack)

    if kind == "radial":
        return solve_grid_timeseries_bfs_fast(
            grid,
            slack_bus_id=slack_bus_id,
            timestep_hours=timestep_hours,
            pf_load_default=pf_load_default,
            load_lagging_default=load_lagging_default,
            pf_pv_default=pf_pv_default,
            pv_lagging_default=pv_lagging_default,
            max_iter=bfs_max_iter,
            tol=bfs_tol,
        )


    if kind == "weakly_meshed":
        g_tree = make_grid_with_lines_removed(grid, drop_line_ids or [])
        res0 = solve_grid_timeseries_bfs(
            g_tree,
            slack_bus_id=slack_bus_id,
            timestep_hours=timestep_hours,
            pf_load_default=pf_load_default,
            load_lagging_default=load_lagging_default,
            pf_pv_default=pf_pv_default,
            pv_lagging_default=pv_lagging_default,
            max_iter=bfs_max_iter,
            tol=bfs_tol,
        )
        return solve_grid_timeseries_newton_raphson(
            grid,
            slack_bus_id=slack_bus_id,
            timestep_hours=timestep_hours,
            pf_load_default=pf_load_default,
            load_lagging_default=load_lagging_default,
            pf_pv_default=pf_pv_default,
            pv_lagging_default=pv_lagging_default,
            S_base_va=nr_S_base_va,
            tol_pu=nr_tol_pu,
            max_iter=nr_max_iter,
            rebuild_J_every=nr_rebuild_J_every,
            Vm0=res0.Vm_pu[0, :] if res0.Vm_pu is not None else None,
            Va0=res0.Va_rad[0, :] if res0.Va_rad is not None else None,
        )

    return solve_grid_timeseries_newton_raphson(
        grid,
        slack_bus_id=slack_bus_id,
        timestep_hours=timestep_hours,
        pf_load_default=pf_load_default,
        load_lagging_default=load_lagging_default,
        pf_pv_default=pf_pv_default,
        pv_lagging_default=pv_lagging_default,
        S_base_va=nr_S_base_va,
        tol_pu=nr_tol_pu,
        max_iter=nr_max_iter,
        rebuild_J_every=nr_rebuild_J_every,
    )

def compute_signed_line_currents(grid, res, *, use_tree=True):
    import numpy as np

    Vmag = np.asarray(res.V, dtype=float)
    Vm_pu = getattr(res, "Vm_pu", None)
    Va = getattr(res, "Va_rad", None)

    slack_bus = int(getattr(res, "slack_bus", getattr(getattr(grid, "kind", None), "slack_bus", 0)))
    node_index = res.node_index
    line_index = res.line_index

    if Vm_pu is not None and Va is not None:
        Vm_pu = np.asarray(Vm_pu, dtype=float)
        Va = np.asarray(Va, dtype=float)
        Vslack_V = float(np.nanmax(Vmag[:, node_index.get(slack_bus, 0)]))
        V_base = Vslack_V if np.isfinite(Vslack_V) and Vslack_V > 0 else 1.0
        Vc = (Vm_pu * np.exp(1j * Va)) * V_base
    else:
        Vc = None

    node_ids = [int(n.bus_id) for n in grid.nodes]
    tree_line_ids = getattr(getattr(grid, "kind", None), "tree_line_ids", None)

    lines_use = list(grid.lines)
    if use_tree and tree_line_ids:
        keep = set(int(x) for x in tree_line_ids)
        lines_use = [l for l in grid.lines if int(l.line_id) in keep]

    adj = {nid: [] for nid in node_ids}
    for l in lines_use:
        u = int(l.from_bus); v = int(l.to_bus)
        if u in adj and v in adj:
            w = float(getattr(l, "length_m", 0.0) or 0.0)
            if not np.isfinite(w) or w < 0: w = 0.0
            adj[u].append((v, w))
            adj[v].append((u, w))

    dist = {nid: np.nan for nid in node_ids}
    if slack_bus not in dist and node_ids:
        slack_bus = node_ids[0]
    dist[slack_bus] = 0.0
    q = [slack_bus]
    seen = {slack_bus}
    while q:
        a = q.pop(0)
        da = float(dist.get(a, 0.0))
        for b, w in adj.get(a, []):
            if b in seen:
                continue
            seen.add(b)
            dist[b] = da + float(w)
            q.append(b)

    T = int(Vmag.shape[0])
    L = len(grid.lines)
    I_signed = np.full((T, L), np.nan, dtype=float)

    for l in grid.lines:
        lid = int(l.line_id)
        k = line_index.get(lid, None)
        if k is None:
            continue

        u_bus = int(l.from_bus)
        v_bus = int(l.to_bus)
        iu = node_index.get(u_bus, None)
        iv = node_index.get(v_bus, None)
        if iu is None or iv is None:
            continue

        rpkm = float(getattr(l, "r_ohm_per_km", 0.0) or 0.0)
        xpkm = float(getattr(l, "x_ohm_per_km", 0.0) or 0.0)
        z = (rpkm + 1j * xpkm) * (float(getattr(l, "length_m", 0.0)) / 1000.0)
        if abs(z) <= 0:
            continue

        if Vc is None:
            continue

        Vu = Vc[:, iu]
        Vv = Vc[:, iv]
        Iuv = (Vu - Vv) / z
        Suv = Vu * np.conj(Iuv)
        Puv = Suv.real

        du = dist.get(u_bus, np.nan)
        dv = dist.get(v_bus, np.nan)

        if np.isfinite(du) and np.isfinite(dv) and du <= dv:
            sgn = np.sign(Puv)
        elif np.isfinite(du) and np.isfinite(dv):
            sgn = -np.sign(Puv)
        else:
            sgn = np.sign(Puv)

        I_signed[:, k] = np.abs(Iuv) * sgn

    return I_signed



def plot_interactive(grid, *, t0=0, results_attr="results"):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    from matplotlib.collections import LineCollection
    import warnings
    import math

    res = getattr(grid, results_attr, None)
    if res is None:
        warnings.warn(f"grid has no '{results_attr}' attribute")
        return
    def _nice_round_up(x):
        x = float(x)
        if not np.isfinite(x) or x <= 0:
            return 1.0
        p = 10 ** math.floor(math.log10(x))
        return math.ceil(x / p) * p
    V = getattr(res, "V", None)
    P_inj = getattr(res, "P_inj", None)
    Q_inj = getattr(res, "Q_inj", None)
    I_line = getattr(res, "I_line", None)
    P_loss_line = getattr(res, "P_loss_line", None)

    if V is None or P_inj is None or Q_inj is None or I_line is None or P_loss_line is None:
        warnings.warn("missing arrays on results")
        return

    V = np.asarray(V, dtype=float)
    P_inj = np.asarray(P_inj, dtype=float)
    Q_inj = np.asarray(Q_inj, dtype=float)
    I_line = np.asarray(I_line, dtype=float)
    P_loss_line = np.asarray(P_loss_line, dtype=float)
    P_all = P_inj
    Q_all = Q_inj
    Pl_all = P_loss_line
    
    P_pos_max = float(np.nanmax(np.sum(np.maximum(P_all, 0.0), axis=1)))
    P_neg_max = float(np.nanmax(np.sum(-np.minimum(P_all, 0.0), axis=1) + np.nansum(Pl_all, axis=1)))
    
    Q_pos_max = float(np.nanmax(np.sum(np.maximum(Q_all, 0.0), axis=1)))
    Q_neg_max = float(np.nanmax(np.sum(-np.minimum(Q_all, 0.0), axis=1)))
    
    need = max(P_pos_max, P_neg_max, Q_pos_max, Q_neg_max, 1e-6) / 1000.0
    PQLIM = _nice_round_up(1.2 * need)
    T, N = V.shape
    t0 = int(np.clip(int(t0), 0, max(T - 1, 0)))

    node_ids = [int(n.bus_id) for n in grid.nodes]
    line_ids = [int(l.line_id) for l in grid.lines]

    node_index = getattr(res, "node_index", {nid: i for i, nid in enumerate(node_ids)})
    line_index = getattr(res, "line_index", {lid: i for i, lid in enumerate(line_ids)})

    slack_bus = int(
        getattr(
            res,
            "slack_bus",
            getattr(getattr(grid, "kind", None), "slack_bus", node_ids[0] if node_ids else 0),
        )
    )



    def _compute_distances_no_heap():
        node_set = set(node_ids)
        adj = {nid: [] for nid in node_ids}
        for l in grid.lines:
            u = int(l.from_bus)
            v = int(l.to_bus)
            if u not in node_set or v not in node_set:
                continue
            w = float(getattr(l, "length_m", 0.0) or 0.0)
            if not np.isfinite(w) or w < 0:
                w = 0.0
            adj[u].append((v, w))
            adj[v].append((u, w))

        dist = {nid: np.inf for nid in node_ids}
        s = slack_bus
        if s not in dist and node_ids:
            s = node_ids[0]
        dist[s] = 0.0

        unvisited = set(node_ids)
        while unvisited:
            u = min(unvisited, key=lambda n: dist[n])
            if not np.isfinite(dist[u]):
                break
            unvisited.remove(u)
            du = dist[u]
            for v, w in adj.get(u, []):
                if v in unvisited:
                    nd = du + w
                    if nd < dist[v]:
                        dist[v] = nd

        for k in dist:
            if not np.isfinite(dist[k]):
                dist[k] = np.nan
        return dist

    dist = _compute_distances_no_heap()

    node_id_arr = np.asarray(node_ids, dtype=int)
    node_dist_arr = np.asarray([float(dist.get(int(nid), np.nan)) for nid in node_id_arr], dtype=float)
    node_col_arr = np.asarray([int(node_index.get(int(nid), -1)) for nid in node_id_arr], dtype=int)

    node_type_arr = np.asarray([str(n.node_type).lower() for n in grid.nodes], dtype=object)
    visible_mask = (node_type_arr != "connection") & np.isfinite(node_dist_arr) & (node_col_arr >= 0)
    allfinite_mask = np.isfinite(node_dist_arr) & (node_col_arr >= 0)

    x_min = float(np.nanmin(node_dist_arr)) if np.any(np.isfinite(node_dist_arr)) else 0.0
    x_max = float(np.nanmax(node_dist_arr)) if np.any(np.isfinite(node_dist_arr)) else 1.0

    Imax_line = np.nanmax(I_line, axis=0) if I_line.ndim == 2 and I_line.shape[1] > 0 else np.array([])
    Imax_line = np.where(np.isfinite(Imax_line) & (Imax_line > 0), Imax_line, 1.0)

    tree_line_ids = getattr(getattr(grid, "kind", None), "tree_line_ids", None)
    if tree_line_ids is not None and len(tree_line_ids) > 0:
        keep = set(int(x) for x in tree_line_ids)
        lines_for_I = [l for l in grid.lines if int(l.line_id) in keep]
    else:
        lines_for_I = list(grid.lines)

    segs_x = []
    seg_line_cols = []
    for l in lines_for_I:
        u = int(l.from_bus)
        v = int(l.to_bus)
        du = dist.get(u, np.nan)
        dv = dist.get(v, np.nan)
        if not (np.isfinite(du) and np.isfinite(dv)):
            continue
        if abs(dv - du) < 1e-12:
            continue
        a, b = (du, dv) if du < dv else (dv, du)
        k = line_index.get(int(l.line_id), None)
        if k is None:
            continue
        segs_x.append((a, b))
        seg_line_cols.append(int(k))

    seg_line_cols = np.asarray(seg_line_cols, dtype=int)
    segs_x = list(segs_x)

    edge_u_cols = []
    edge_v_cols = []
    edge_du = []
    edge_dv = []
    for l in grid.lines:
        u = int(l.from_bus)
        v = int(l.to_bus)
        du = dist.get(u, np.nan)
        dv = dist.get(v, np.nan)
        iu = node_index.get(u, None)
        iv = node_index.get(v, None)
        if iu is None or iv is None:
            continue
        if not (np.isfinite(du) and np.isfinite(dv)):
            continue
        edge_u_cols.append(int(iu))
        edge_v_cols.append(int(iv))
        edge_du.append(float(du))
        edge_dv.append(float(dv))

    edge_u_cols = np.asarray(edge_u_cols, dtype=int)
    edge_v_cols = np.asarray(edge_v_cols, dtype=int)
    edge_du = np.asarray(edge_du, dtype=float)
    edge_dv = np.asarray(edge_dv, dtype=float)

    def _voltage_colors(v):
        v = np.asarray(v, dtype=float)
        c = np.full(v.shape, "black", dtype=object)
        c[(v < 360.0) | (v > 440.0)] = "red"
        c[((v >= 360.0) & (v < 380.0)) | ((v > 420.0) & (v <= 440.0))] = "orange"
        return c

    def _current_colors(i_pu):
        i_pu = np.asarray(i_pu, dtype=float)
        c = np.full(i_pu.shape, "black", dtype=object)
        c[i_pu > 1.0] = "red"
        c[(i_pu > 0.9) & (i_pu <= 1.0)] = "orange"
        return c

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.1, 0.9], hspace=0.35, wspace=0.25)

    axV = fig.add_subplot(gs[0, 0])
    axI = fig.add_subplot(gs[0, 1])
    P_pos_max = float(np.nanmax(np.sum(np.maximum(P_inj, 0.0), axis=1))) / 1000.0
    P_abs_max = float(np.nanmax(np.sum(-np.minimum(P_inj, 0.0), axis=1))) / 1000.0
    P_loss_max = float(np.nanmax(np.nansum(P_loss_line, axis=1))) / 1000.0
    
    Q_pos_max = float(np.nanmax(np.sum(np.maximum(Q_inj, 0.0), axis=1))) / 1000.0
    Q_abs_max = float(np.nanmax(np.sum(-np.minimum(Q_inj, 0.0), axis=1))) / 1000.0
    
    pos_need = max(P_pos_max, Q_pos_max, 1e-6)
    neg_need = max(P_abs_max + P_loss_max, Q_abs_max, 1e-6)
    
    ylim_pos = _nice_round_up(1.2 * pos_need)
    ylim_neg = -_nice_round_up(1.2 * neg_need)
    axPQ = fig.add_subplot(gs[1, :])
    axTxt = fig.add_subplot(gs[2, :])
    axTxt.axis("off")

    ax_slider = fig.add_axes([0.10, 0.03, 0.80, 0.04])
    slider = Slider(ax_slider, "t", 0, max(T - 1, 0), valinit=t0, valstep=1)
    fig._pf_slider = slider

    edge_lc = LineCollection([], linewidths=0.3, colors="black")
    axV.add_collection(edge_lc)

    V0_all = V[t0, :]
    y_nodes0 = V0_all[node_col_arr]
    sizes = np.zeros_like(y_nodes0, dtype=float)
    sizes[visible_mask] = 22.0

    colors0 = _voltage_colors(y_nodes0)
    colors0[~visible_mask] = "none"

    scV = axV.scatter(node_dist_arr, y_nodes0, s=sizes, c=colors0, edgecolors="none")

    axV.set_xlim(x_min, x_max*1.1)
    axV.set_ylim(340.0, 460.0)
    axV.set_xlabel("Distance (m)")
    axV.set_ylabel("Voltage (V)")
    axV.set_title("")

    lcI = LineCollection([], linewidths=2.0)
    axI.add_collection(lcI)
    axI.set_xlim(x_min, x_max*1.1)
    axI.set_ylim(0.0, 2.0)
    axI.set_xlabel("Distance  (m)")
    axI.set_ylabel("I / Imax_line (pu)")
    axI.set_title("")

    xcats = np.array([0.0, 1.0])
    width = 0.55

    bP_in = axPQ.bar([xcats[0]], [0.0], width=width, color="green", alpha=1.0)
    bP_abs = axPQ.bar([xcats[0]], [0.0], width=width, color="gray", alpha=1.0)
    bP_loss = axPQ.bar([xcats[0]], [0.0], width=width, color="red", alpha=1.0)

    bQ_in = axPQ.bar([xcats[1]], [0.0], width=width, color="green", alpha=0.6)
    bQ_abs = axPQ.bar([xcats[1]], [0.0], width=width, color="gray", alpha=0.6)

    axPQ.set_xticks([0.0, 1.0], ["P (kW)", "Q (kVAr)"])
    axPQ.axhline(0.0, linewidth=1.2)
    axPQ.set_ylabel("Power")
    axPQ.set_ylim(ylim_neg, ylim_pos)
    txt = axTxt.text(0.01, 0.95, "", va="top", family="monospace")

    def _set_bar(bar_container, val):
        for rect in bar_container:
            rect.set_height(float(val))

    def _update_voltage_edges(t):
        Vt = V[t, :]
        Vu = Vt[edge_u_cols]
        Vv = Vt[edge_v_cols]
        segs = [((float(du), float(vu)), (float(dv), float(vv))) for du, dv, vu, vv in zip(edge_du, edge_dv, Vu, Vv)]
        edge_lc.set_segments(segs)

    def _update_segments_I(t):
        if len(segs_x) == 0:
            lcI.set_segments([])
            return
        Ipu = (I_line[t, seg_line_cols] / Imax_line[seg_line_cols]).astype(float)
        cols = _current_colors(Ipu)
        segs = [((a, float(i)), (b, float(i))) for (a, b), i in zip(segs_x, Ipu)]
        lcI.set_segments(segs)
        lcI.set_color(list(cols))

    def _update_power_axis(t):
        P = P_inj[t, :]
        Q = Q_inj[t, :]
        Pl = P_loss_line[t, :]

        P_in = float(np.sum(P[P > 0.0])) / 1000.0 if P.size else 0.0
        P_abs = float(-np.sum(P[P < 0.0])) / 1000.0 if P.size else 0.0
        P_loss = float(np.nansum(Pl)) / 1000.0 if Pl.size else 0.0

        Q_in = float(np.sum(Q[Q > 0.0])) / 1000.0 if Q.size else 0.0
        Q_abs = float(-np.sum(Q[Q < 0.0])) / 1000.0 if Q.size else 0.0

        P_abs_neg = -P_abs
        P_loss_neg = -P_loss
        Q_abs_neg = -Q_abs

        _set_bar(bP_in, P_in)
        _set_bar(bP_abs, P_abs_neg)
        _set_bar(bP_loss, P_loss_neg)
        _set_bar(bQ_in, Q_in)
        _set_bar(bQ_abs, Q_abs_neg)

        need_pos = max(P_in, Q_in, 0.0)
        need_neg = max(abs(P_abs_neg + P_loss_neg), abs(Q_abs_neg), 0.0)


        Vt = V[t, :]
        Vmin = float(np.nanmin(Vt)) if Vt.size else np.nan
        Vmax = float(np.nanmax(Vt)) if Vt.size else np.nan
        Ipu_max = float(np.nanmax(I_line[t, :] / Imax_line)) if I_line.shape[1] else np.nan

        txt.set_text(
            "\n".join(
                [
                    f"t = {t}",
                    f"slack_bus = {slack_bus}",
                    f"Vmin = {Vmin:.3f} V   Vmax = {Vmax:.3f} V",
                    f"max(I/Imax_line) = {Ipu_max:.3f}",
                    f"P_inj = {P_in:.3f} kW   P_abs = {P_abs:.3f} kW   P_loss = {P_loss:.3f} kW",
                    f"Q_inj = {Q_in:.3f} kVAr  Q_abs = {Q_abs:.3f} kVAr",
                ]
            )
        )

    def _update(t):
        t = int(t)
        Vt = V[t, :]
        y_nodes = Vt[node_col_arr]

        cols = _voltage_colors(y_nodes)
        cols[~visible_mask] = "none"

        scV.set_offsets(np.c_[node_dist_arr, y_nodes])
        scV.set_sizes(np.where(visible_mask, 22.0, 0.0))
        scV.set_color(cols)

        _update_voltage_edges(t)
        _update_segments_I(t)
        _update_power_axis(t)

        fig.canvas.draw_idle()

    def _update_safe(val):
        try:
            _update(int(round(val)))
        except Exception as e:
            print("plot_interactive update error:", repr(e))

    fig._pf_update = _update_safe
    slider.on_changed(_update_safe)

    _update(t0)
    plt.show(block=True)
    return fig


from pathlib import Path


def generate_grid_fault_report(
    gs,
    min_voltage,
    max_voltage,
    max_supply_apparent_power,
    out_path="grid_fault_report.txt",
):
    out_path = Path(out_path)

    def to_1d(arr):
        return np.asarray(arr).reshape(-1)

    with out_path.open("w", encoding="utf-8") as f:
        for i, grid in enumerate(gs.grids):
            f.write(f"Grid {i}\n\n")
            f.write("Faults:\n")
            V_max = 0,
            V_min = 1000
            f.write("Overvoltage:\n")
            for ni, node in enumerate(grid.nodes):
                if node.node_type != "building":
                    continue
                V = to_1d(node.V)
                steps = np.where(V > max_voltage)[0]
                if steps.size > 0:
                    f.write(f"node{ni} at {steps.tolist()}\n")
                if max(V)>V_max:
                    V_max = max(V)
            
            f.write(f"maximum voltage is {V_max}\n")
            f.write("\nUndervoltage:\n")
            for ni, node in enumerate(grid.nodes):
                if node.node_type != "building":
                    continue
                V = to_1d(node.V)
                steps = np.where(V < min_voltage)[0]
                if steps.size > 0:
                    f.write(f"node{ni} at {steps.tolist()}\n")
                    if min(V)<V_min:
                        V_min = min(V)
            f.write(f"minimum voltage is {V_min}\n")
            f.write("\nOvercurrent:\n")
            for li, line in enumerate(grid.lines):
                I = to_1d(line.I)
                steps = np.where(I > line.max_current_a)[0]
                if steps.size > 0:
                    f.write(f"line{li} at {steps.tolist()}\n")
                
            f.write("\nOverflow:\n")
            S = to_1d(grid.results.S_supply)
            steps = np.where(S > max_supply_apparent_power)[0]
            if steps.size > 0:
                f.write(f"at {steps.tolist()}\n")

            f.write("\nBackflow:\n")
            P = to_1d(grid.results.P_supply)
            steps = np.where(P < 0)[0]
            if steps.size > 0:
                f.write(f"at {steps.tolist()}\n")
            

            f.write("\n" + "-" * 40 + "\n\n")

    return out_path



