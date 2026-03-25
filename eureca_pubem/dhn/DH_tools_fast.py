from __future__ import annotations

import numpy as np

try:
    from numba import njit
    _HAS_NUMBA = True
except Exception:
    _HAS_NUMBA = False


def _as1d(x) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        return a.reshape(1)
    return a.reshape(-1)


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


def _build_m_ts(system, t_idx: np.ndarray) -> np.ndarray:
    lines = system.lines
    E = len(lines)
    T = int(t_idx.shape[0])
    m_ts = np.empty((T, E), dtype=float)
    for e, l in enumerate(lines):
        m = getattr(l, "m_init", None)
        if m is None:
            m = getattr(l, "m", None)
        if m is None:
            raise ValueError(f"Line {getattr(l,'line_id',e)} missing both m_init and m.")
        m = _as1d(m)
        if m.shape[0] == 1:
            m_ts[:, e] = float(m[0])
        else:
            m_ts[:, e] = m[t_idx]
    return m_ts


def _build_demand_ts(system, t_idx: np.ndarray) -> np.ndarray:
    nodes = system.nodes
    N = len(nodes)
    T = int(t_idx.shape[0])
    dem = np.zeros((T, N), dtype=float)
    for i, n in enumerate(nodes):
        d = getattr(n, "demand", None)
        if d is None:
            continue
        d = _as1d(d)
        if d.shape[0] == 1:
            dem[:, i] = float(d[0])
        else:
            dem[:, i] = d[t_idx]
    return dem


def _node_type_vec(system) -> np.ndarray:
    nodes = system.nodes
    out = np.zeros(len(nodes), dtype=np.int64)
    for i, n in enumerate(nodes):
        tp = getattr(n, "node_type", None)
        if tp == "consumer":
            out[i] = 1
        elif tp == "supply":
            out[i] = 2
        else:
            out[i] = 0
    return out


def _supply_pos(system) -> int:
    for i, n in enumerate(system.nodes):
        if getattr(n, "node_type", None) == "supply":
            return i
    raise ValueError("No supply node.")


def _eff_vec(system, T: int, t_idx: np.ndarray, eff_default=0.99, attr_eff="node_dhn_efficiency") -> np.ndarray:
    nodes = system.nodes
    eff = np.full((T, len(nodes)), float(eff_default), dtype=float)
    for i, n in enumerate(nodes):
        e = getattr(n, attr_eff, None)
        if e is None:
            continue
        e = _as1d(e)
        if e.shape[0] == 1:
            eff[:, i] = float(e[0])
        else:
            eff[:, i] = e[t_idx]
    return eff


if _HAS_NUMBA:
    @njit(cache=False)
    def _pressures_timeseries_nb(
        m_ts,
        L,
        D,
        eps_m,
        parent,
        parent_edge,
        order,
        sign,
        datum,
        rho,
        mu_supply,
        mu_return,
        return_pressure_default,
        required_dp_default,
    ):
        T, E = m_ts.shape
        N = parent.shape[0]

        dp_sup = np.empty((T, E), dtype=np.float64)
        dp_ret = np.empty((T, E), dtype=np.float64)
        Ps = np.empty((T, N), dtype=np.float64)
        Pr = np.empty((T, N), dtype=np.float64)

        A = np.pi * D * D / 4.0

        i = 0
        for t in range(T):
            for e in range(E):
                i = i+1
                m = m_ts[t, e]
                v = m / (rho * A[e])
                vabs = abs(v)

                Re = rho * vabs * D[e] / mu_supply
                if Re < 1e-12:
                    Re = 1e-12
                if Re < 2300.0:
                    f = 64.0 / Re
                else:
                    f = 0.25 / (np.log10(eps_m[e] / (3.7 * D[e]) + 5.74 / (Re ** 0.9)) ** 2)
                dp_sup[t, e] = f * (L[e] / D[e]) * 0.5 * rho * v * vabs

                Re2 = rho * vabs * D[e] / mu_return
                if Re2 < 1e-12:
                    Re2 = 1e-12
                if Re2 < 2300.0:
                    f2 = 64.0 / Re2
                else:
                    f2 = 0.25 / (np.log10(eps_m[e] / (3.7 * D[e]) + 5.74 / (Re2 ** 0.9)) ** 2)
                dp_ret[t, e] = f2 * (L[e] / D[e]) * 0.5 * rho * v * vabs

            # print(i, T, E)
            for i in range(N):
                Ps[t, i] = np.nan
                Pr[t, i] = np.nan

            Ps[t, datum] = 0.0
            for k in range(order.shape[0]):
                vpos = order[k]
                if vpos == datum:
                    continue
                ppos = parent[vpos]
                eidx = parent_edge[vpos]
                Ps[t, vpos] = Ps[t, ppos] + sign[vpos] * dp_sup[t, eidx]

            Pr[t, datum] = return_pressure_default
            for k in range(order.shape[0]):
                vpos = order[k]
                if vpos == datum:
                    continue
                ppos = parent[vpos]
                eidx = parent_edge[vpos]
                Pr[t, vpos] = Pr[t, ppos] - sign[vpos] * dp_ret[t, eidx]

            min_delta = 1e99
            for i in range(N):
                dlt = Ps[t, i] - Pr[t, i]
                if dlt < min_delta:
                    min_delta = dlt

            shift = required_dp_default - min_delta
            for i in range(N):
                Ps[t, i] = Ps[t, i] + shift

        return Ps, Pr, dp_sup, dp_ret


    @njit(cache=False)
    def _thermal_timeseries_tree_nb(
        m_ts,
        demand_kw_ts,
        eff_ts,
        UA,
        earth_T_ts,
        parent,
        parent_edge,
        children_ptr,
        children_idx,
        order,
        sign,
        node_type,
        supply_pos,
        Ts_supply_set,
        cp,
    ):
        T, E = m_ts.shape
        N = parent.shape[0]

        Ts_node = np.empty((T, N), dtype=np.float64)
        Tr_node = np.empty((T, N), dtype=np.float64)

        Ts_start = np.empty((T, E), dtype=np.float64)
        Ts_end = np.empty((T, E), dtype=np.float64)
        Tr_start = np.empty((T, E), dtype=np.float64)
        Tr_end = np.empty((T, E), dtype=np.float64)

        for t in range(T):
            for i in range(N):
                Ts_node[t, i] = np.nan
                Tr_node[t, i] = np.nan

            Ts_node[t, supply_pos] = Ts_supply_set

            for k in range(order.shape[0]):
                u = order[k]
                tu = Ts_node[t, u]
                if not np.isfinite(tu):
                    continue

                c0 = children_ptr[u]
                c1 = children_ptr[u + 1]
                for jj in range(c0, c1):
                    v = children_idx[jj]
                    e = parent_edge[v]

                    mdot = abs(m_ts[t, e])
                    if mdot < 1e-12:
                        mdot = 1e-12

                    Tin = tu
                    Te = earth_T_ts[t]

                    Ts_start[t, e] = Tin
                    alpha = np.exp(-UA[e] / (mdot * cp))
                    Tout = Te + (Tin - Te) * alpha
                    Ts_end[t, e] = Tout

                    Ts_node[t, v] = Tout

            for i in range(N):
                if node_type[i] == 1:
                    P_w = demand_kw_ts[t, i] * 1000.0
                    mdot_n = 0.0
                    c0 = children_ptr[i]
                    c1 = children_ptr[i + 1]
                    if c0 == c1:
                        e = parent_edge[i]
                        mdot_n = abs(m_ts[t, e])
                    else:
                        s = 0.0
                        for jj in range(c0, c1):
                            v = children_idx[jj]
                            e = parent_edge[v]
                            s += abs(m_ts[t, e])
                        mdot_n = s
                    if mdot_n < 1e-12:
                        mdot_n = 1e-12
                    dT = P_w / (mdot_n * cp * eff_ts[t, i])
                    Tr_node[t, i] = Ts_node[t, i] - dT

            for kk in range(order.shape[0] - 1, -1, -1):
                u = order[kk]
                if node_type[u] == 1:
                    continue

                c0 = children_ptr[u]
                c1 = children_ptr[u + 1]
                if c0 == c1:
                    continue

                msum = 0.0
                Tmix = 0.0
                for jj in range(c0, c1):
                    v = children_idx[jj]
                    e = parent_edge[v]

                    mdot = abs(m_ts[t, e])
                    if mdot < 1e-12:
                        mdot = 1e-12

                    Tv = Tr_node[t, v]
                    if not np.isfinite(Tv):
                        continue

                    Te = earth_T_ts[t]
                    alpha = np.exp(-UA[e] / (mdot * cp))
                    Tin = Te + (Tv - Te) * alpha

                    Tr_end[t, e] = Tv
                    Tr_start[t, e] = Tin

                    msum += mdot
                    Tmix += mdot * Tin

                if msum > 0.0:
                    Tr_node[t, u] = Tmix / msum

            for i in range(N):
                if i == supply_pos:
                    if not np.isfinite(Tr_node[t, i]):
                        Tr_node[t, i] = Ts_node[t, i] - 20.0

        return Ts_node, Tr_node, Ts_start, Ts_end, Tr_start, Tr_end


def _children_csr(cache_children) -> tuple[np.ndarray, np.ndarray]:
    N = len(cache_children)
    ptr = np.zeros(N + 1, dtype=np.int64)
    nnz = 0
    for i in range(N):
        nnz += len(cache_children[i])
        ptr[i + 1] = nnz
    idx = np.zeros(nnz, dtype=np.int64)
    k = 0
    for i in range(N):
        for v in cache_children[i]:
            idx[k] = int(v)
            k += 1
    return ptr, idx


def solve_timeseries_fast_tree(
    system,
    time_frame,
    *,
    return_pressure_default=150_000.0,
    required_pressure_difference_default=100_000.0,
    supply_temperature=80.0,

    min_supply_temperature=None,
    supply_temperature_margin=0.5,
    max_supply_temperature=50000.0,
    supply_temperature_tol=0.05,

    cp=4180.0,
    rho=1000.0,
    mu_supply=0.0003372,
    mu_return=0.0008891,
    eff_default=0.99,
    attr_eff="node_dhn_efficiency",
):
    if not _HAS_NUMBA:
        raise RuntimeError("Numba not available.")

    # ------------------------------------------------------------
    # Geometry & topology
    # ------------------------------------------------------------
    cache = system.pressure_cache

    parent = np.asarray(cache["parent"], dtype=np.int64)
    parent_edge = np.asarray(cache["parent_edge"], dtype=np.int64)
    order = np.asarray(cache["order"], dtype=np.int64)
    sign = np.asarray(cache["sign"], dtype=np.float64)
    datum = int(cache["datum_pos"])

    L = np.asarray(cache["L"], dtype=np.float64)
    D = np.asarray(cache["D"], dtype=np.float64)
    eps = np.asarray(cache["eps"], dtype=np.float64) / 1000.0

    children_ptr, children_idx = _children_csr(cache["children"])

    # ------------------------------------------------------------
    # Time handling
    # ------------------------------------------------------------
    t_idx = np.asarray(list(time_frame), dtype=np.int64)
    T = t_idx.size
    if T == 0:
        return

    # ------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------
    m_ts = _build_m_ts(system, t_idx)
    dem_ts = _build_demand_ts(system, t_idx)
    eff_ts = _eff_vec(system, T, t_idx, eff_default, attr_eff)

    lines = system.lines
    UA = np.array([float(l.UA) for l in lines])

    nodes = system.nodes
    N = len(nodes)
    node_type = _node_type_vec(system)
    supply_pos = _supply_pos(system)

    # ------------------------------------------------------------
    # Earth temperature
    # ------------------------------------------------------------
    earth_T = cache.get("earth_temperature", None)
    if earth_T is None:
        earth_T = getattr(system, "ground_temperature_array", None)
    if earth_T is None:
        raise ValueError("No earth temperature array found.")

    earth_T = _as1d(earth_T)
    earth_T_ts = earth_T[t_idx] if earth_T.size > 1 else np.full(T, earth_T[0])

    # ------------------------------------------------------------
    # PRESSURE SOLVE (ONCE)
    # ------------------------------------------------------------
    Ps, Pr, dp_sup, dp_ret = _pressures_timeseries_nb(
        m_ts,
        L,
        D,
        eps,
        parent,
        parent_edge,
        order,
        sign,
        datum,
        float(rho),
        float(mu_supply),
        float(mu_return),
        float(return_pressure_default),
        float(required_pressure_difference_default),
    )

    # ------------------------------------------------------------
    # PER-TIMESTEP SUPPLY TEMPERATURE
    # ------------------------------------------------------------
    Ts_supply_ts = np.full(T, float(supply_temperature))

    mask_cons = node_type == 1

    for ti in range(T):

        Ts0 = Ts_supply_ts[ti]

        Ts_node_1, _, _, _, _, _ = _thermal_timeseries_tree_nb(
            m_ts[ti:ti+1],
            dem_ts[ti:ti+1],
            eff_ts[ti:ti+1],
            UA,
            earth_T_ts[ti:ti+1],
            parent,
            parent_edge,
            children_ptr,
            children_idx,
            order,
            sign,
            node_type,
            int(supply_pos),
            Ts0,
            float(cp),
        )

        if min_supply_temperature is None:
            continue

        Ts_cons = Ts_node_1[0, mask_cons]
        Tmin = np.nanmin(Ts_cons)

        if Tmin >= min_supply_temperature - supply_temperature_tol:
            continue

        Te = earth_T_ts[ti]

        Gamma = (Tmin - Te) / (Ts0 - Te)
        if Gamma <= 0.0 or not np.isfinite(Gamma):
            raise RuntimeError(
                f"Invalid attenuation at timestep {ti}: Γ={Gamma}"
            )

        Ts_req = Te + (min_supply_temperature - Te) / Gamma
        Ts_req += supply_temperature_margin

        if Ts_req > max_supply_temperature:
            raise RuntimeError(
                f"Required Ts={Ts_req:.1f} °C exceeds maximum at t={ti}"
            )

        Ts_supply_ts[ti] = Ts_req

    # ------------------------------------------------------------
    # FINAL THERMAL SOLVE (PER TIMESTEP)
    # ------------------------------------------------------------
    Ts_node = np.zeros((T, N))
    Tr_node = np.zeros((T, N))
    Ts_start = np.zeros((T, len(lines)))
    Ts_end = np.zeros((T, len(lines)))
    Tr_start = np.zeros((T, len(lines)))
    Tr_end = np.zeros((T, len(lines)))

    for ti in range(T):
        Ts_n, Tr_n, Ts_s, Ts_e, Tr_s, Tr_e = _thermal_timeseries_tree_nb(
            m_ts[ti:ti+1],
            dem_ts[ti:ti+1],
            eff_ts[ti:ti+1],
            UA,
            earth_T_ts[ti:ti+1],
            parent,
            parent_edge,
            children_ptr,
            children_idx,
            order,
            sign,
            node_type,
            int(supply_pos),
            Ts_supply_ts[ti],
            float(cp),
        )

        Ts_node[ti] = Ts_n[0]
        Tr_node[ti] = Tr_n[0]
        Ts_start[ti] = Ts_s[0]
        Ts_end[ti] = Ts_e[0]
        Tr_start[ti] = Tr_s[0]
        Tr_end[ti] = Tr_e[0]

    # ------------------------------------------------------------
    # WRITE BACK TO OBJECTS
    # ------------------------------------------------------------
    Tfull = int(cache["T"])

    for n in nodes:
        _ensure_arr(n, "supply_pressure", Tfull)
        _ensure_arr(n, "return_pressure", Tfull)
        _ensure_arr(n, "supply_temperature", Tfull)
        _ensure_arr(n, "return_temperature", Tfull)

    for l in lines:
        _ensure_arr(l, "supply_pressure_drop", Tfull)
        _ensure_arr(l, "return_pressure_drop", Tfull)
        _ensure_arr(l, "supply_start_temperature", Tfull)
        _ensure_arr(l, "supply_end_temperature", Tfull)
        _ensure_arr(l, "return_start_temperature", Tfull)
        _ensure_arr(l, "return_end_temperature", Tfull)

    for ti, t in enumerate(t_idx):
        for i, n in enumerate(nodes):
            n.supply_pressure[t] = Ps[ti, i]
            n.return_pressure[t] = Pr[ti, i]
            n.supply_temperature[t] = Ts_node[ti, i]
            n.return_temperature[t] = Tr_node[ti, i]

        for e, l in enumerate(lines):
            l.supply_pressure_drop[t] = dp_sup[ti, e]
            l.return_pressure_drop[t] = dp_ret[ti, e]
            l.supply_start_temperature[t] = Ts_start[ti, e]
            l.supply_end_temperature[t] = Ts_end[ti, e]
            l.return_start_temperature[t] = Tr_start[ti, e]
            l.return_end_temperature[t] = Tr_end[ti, e]