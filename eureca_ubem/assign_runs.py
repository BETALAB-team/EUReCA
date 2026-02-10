from __future__ import annotations

from typing import Any, Dict, List, Optional
import numpy as np
import copy

_NO_CHANGE = "__NO_CHANGE__"


def _scenario_key_from_values(values_in_order: List[Any]) -> str:
    """Builds scenario key, skipping NO_CHANGE values; spaces->'_' only in key."""
    parts = [str(v).replace(" ", "_") for v in values_in_order if v != _NO_CHANGE]
    return "scenario_" + ("base" if not parts else "_".join(parts))


def _normalize_probs(p: np.ndarray, dim: str, eps: float = 1e-12) -> np.ndarray:
    if np.any(p < -eps):
        raise ValueError(f"Diffusion probs for '{dim}' contain negative values.")
    p = np.clip(p, 0.0, None)
    s = float(p.sum())
    if s < eps:
        raise ValueError(f"Diffusion probs for '{dim}' sum to ~0; cannot sample.")
    return p / s


def generate_heatnode_diffusion_runs(
    scenario_cities: Dict[str, Any],
    scenario_changes: Dict[str, List[Any]],
    scenario_diffusion: Dict[str, List[float]],
    n_runs: int,
    *,
    rng: Optional[np.random.Generator] = None,
    deep_copy_nodes: bool = False,
    diffusion_includes_no_change: Optional[bool] = None,
    require_all_sampled_scenarios_exist: bool = True,
) -> Dict[int, Dict[Any, Any]]:
    """
    Randomly assigns each heat node to a scenario city, matching diffusion marginals per dimension.

    - Supports "dimension not applied" via a NO_CHANGE option.
    - Scenario keys are like: scenario_A_Heat_Pump (spaces -> underscores in key only)
      and if all dims are NO_CHANGE => scenario_base.

    Parameters
    ----------
    scenario_cities:
        scenario_key -> City object (each has .dhn_nodes: dict[node_id -> HeatNode])
        Must include at least scenario_base, plus any combos you want to sample.
    scenario_changes:
        dimension -> list of *real* options (NO_CHANGE not included here)
        Order is used for mapping probs and building scenario keys.
    scenario_diffusion:
        dimension -> list of probs.
        Two styles:
          A) implicit NO_CHANGE: len == len(options) and p_none = 1 - sum(p)
          B) explicit NO_CHANGE: len == len(options)+1 with p[0]=p_none, p[1:]=options probs
    diffusion_includes_no_change:
        If None: auto-detect per dimension based on length.
        If True/False: enforce that style for all dims.
    require_all_sampled_scenarios_exist:
        If True: error if any sampled scenario key is missing from scenario_cities.
        If False: fallback to scenario_base when missing.

    Returns
    -------
    runs:
        {run_idx: {node_id: HeatNode}} where HeatNode is taken from the assigned scenario city.
    """
    if rng is None:
        rng = np.random.default_rng()

    if not scenario_cities:
        raise ValueError("scenario_cities is empty.")

    if "scenario_base" not in scenario_cities:
        raise KeyError("scenario_cities must include 'scenario_base'.")

    # Decide which dimensions participate = those present in BOTH dicts
    dims = [d for d in scenario_changes.keys() if d in scenario_diffusion]
    if not dims:
        raise ValueError("No overlapping dimensions between scenario_changes and scenario_diffusion.")

    # Validate dhn_nodes consistency across scenarios
    base_city = scenario_cities["scenario_base"]
    if not hasattr(base_city, "dhn_nodes") or not isinstance(base_city.dhn_nodes, dict):
        raise TypeError("scenario_base city has no valid '.dhn_nodes' dict.")

    node_ids = list(base_city.dhn_nodes.keys())
    node_id_set = set(node_ids)

    for skey, city in scenario_cities.items():
        if not hasattr(city, "dhn_nodes") or not isinstance(city.dhn_nodes, dict):
            raise TypeError(f"City for scenario '{skey}' has no valid '.dhn_nodes' dict.")
        if set(city.dhn_nodes.keys()) != node_id_set:
            raise ValueError(f"Scenario '{skey}' has different dhn_nodes ids than scenario_base.")

    n_nodes = len(node_ids)

    # Build per-dimension option lists including NO_CHANGE + probability vectors
    opts_per_dim: Dict[str, List[Any]] = {}
    probs_per_dim: Dict[str, np.ndarray] = {}

    for dim in dims:
        real_opts = list(scenario_changes[dim])
        p_in = np.array(scenario_diffusion[dim], dtype=float)

        # Determine whether p includes NO_CHANGE
        if diffusion_includes_no_change is None:
            if len(p_in) == len(real_opts) + 1:
                includes_none = True
            elif len(p_in) == len(real_opts):
                includes_none = False
            else:
                raise ValueError(
                    f"'{dim}': diffusion length {len(p_in)} must be "
                    f"{len(real_opts)} (implicit NO_CHANGE) or {len(real_opts)+1} (explicit)."
                )
        else:
            includes_none = diffusion_includes_no_change
            expected = len(real_opts) + (1 if includes_none else 0)
            if len(p_in) != expected:
                raise ValueError(
                    f"'{dim}': expected diffusion length {expected} (diffusion_includes_no_change={includes_none}), "
                    f"got {len(p_in)}."
                )

        if includes_none:
            # p_in[0] = p_none, p_in[1:] aligned with real options
            opts = [_NO_CHANGE] + real_opts
            p = p_in
        else:
            # implicit none
            p_none = 1.0 - float(p_in.sum())
            if p_none < -1e-9:
                raise ValueError(
                    f"'{dim}': sum(probs)={float(p_in.sum()):.4f} > 1, cannot infer NO_CHANGE."
                )
            p_none = max(0.0, p_none)
            opts = [_NO_CHANGE] + real_opts
            p = np.concatenate([[p_none], p_in])

        p = _normalize_probs(p, dim)
        opts_per_dim[dim] = opts
        probs_per_dim[dim] = p

    # Monte Carlo runs
    runs: Dict[int, Dict[Any, Any]] = {}

    for run in range(n_runs):
        # sample each dim for each node (values include NO_CHANGE)
        sampled: Dict[str, np.ndarray] = {}
        for dim in dims:
            sampled[dim] = rng.choice(opts_per_dim[dim], size=n_nodes, p=probs_per_dim[dim])

        assigned_nodes: Dict[Any, Any] = {}

        for i, nid in enumerate(node_ids):
            values_in_order = [sampled[dim][i] for dim in dims]
            skey = _scenario_key_from_values(values_in_order)

            if skey not in scenario_cities:
                if require_all_sampled_scenarios_exist:
                    raise KeyError(
                        f"Sampled scenario '{skey}' not found in scenario_cities. "
                        f"Either generate it, or set require_all_sampled_scenarios_exist=False."
                    )
                skey = "scenario_base"

            node_obj = scenario_cities[skey].dhn_nodes[nid]
            if deep_copy_nodes:
                node_obj = copy.deepcopy(node_obj)

            assigned_nodes[nid] = node_obj

        runs[run] = assigned_nodes

    return runs
