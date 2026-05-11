import pandas as pd
import geopandas as gpd
import numpy as np
import json
import itertools
import math
from pathlib import Path
from eureca_pubem import scenario_process as sc


ENV_ORDER = ["none", "shallow", "medium", "deep"]

ENV_ALLOWED = {
    "none": ["none", "shallow", "medium", "deep"],
    "shallow": ["shallow", "medium", "deep"],
    "medium": ["medium", "deep"],
    "deep": ["deep"]
}

HEAT_ALLOWED = {
    "boiler": ["dhn", "hp_le", "hp_me", "hp_he"],
    "dhn": ["dhn", "hp_le", "hp_me", "hp_he"],
    "hp_le": ["dhn", "hp_le", "hp_me", "hp_he"],
    "hp_me": ["dhn", "hp_me", "hp_he"],
    "hp_he": ["hp_he", "dhn"]
}

DHW_ALLOWED = HEAT_ALLOWED.copy()

PV_STEPS = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]

HP_TYPES = {"hp_le", "hp_me", "hp_he"}


def get_heat_key(value):
    if not isinstance(value, str):
        raise ValueError(value)

    if "boiler" in value:
        return "boiler"
    if "dhn" in value:
        return "dhn"
    if "hp_le" in value:
        return "hp_le"
    if "hp_me" in value:
        return "hp_me"
    if "hp_he" in value:
        return "hp_he"

    raise ValueError(value)


def compare_config_with_reference(
    reference_dictionary,
    baseline_gdf_path,
    configuration,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20
):
    interv_dict, building_info = sc.make_dictionary(
        baseline_geojson=baseline_gdf_path,
        city=mycity,
        baseline_scenario=baseline_scenario,
        weatherfile_path=weatherfile_path,
        intervention_dictionary=configuration
    )

    new_dict = build_dict_gen(
        building_info=building_info,
        intervention_dict=interv_dict,
        baseline_gdf_path=baseline_gdf_path,
        configuration=configuration,
        ee_measure_path=ee_measure_path,
        pv_type_path=pv_type_path,
        hp_catalog_path=hp_catalog_path,
        grid_pricing_path=grid_pricing_path,
        dhn_pricing_path=dhn_pricing_path,
        fuels_path=fuels_path,
        spot_price_path=spot_price_path
    )

    annuity_factor = (1 - (1 + r) ** (-T)) / r

    for idx, building in new_dict.items():
        building["costs"]["financial"] = {}

        reference_building = reference_dictionary[idx]

        costs_new_hourly = building["costs"]["operational cost"]["hourly_total"]
        costs_reference_hourly = reference_building["costs"]["operational cost"]["hourly_total"]

        building["costs"]["financial"]["savings_hourly"] = costs_reference_hourly - costs_new_hourly
        building["costs"]["financial"]["savings_yearly"] = np.sum(
            building["costs"]["financial"]["savings_hourly"]
        )

        S = building["costs"]["financial"]["savings_yearly"]
        C = building["costs"]["capital cost"]["total"]

        building["costs"]["financial"]["NPV"] = -C + S * annuity_factor

    return new_dict


def generate_configs_for_building(b0, i=0):
    env_options = ENV_ALLOWED[b0["env"]]

    heat_options = HEAT_ALLOWED[get_heat_key(b0["heat"])].copy()
    if b0["heat"] not in heat_options:
        heat_options.append(b0["heat"])

    dhw_options = DHW_ALLOWED[get_heat_key(b0["dhw"])].copy()
    if b0["dhw"] not in dhw_options:
        dhw_options.append(b0["dhw"])

    pv_min = b0["pv_percentage"]
    pv_options = [p for p in PV_STEPS if p >= pv_min]

    configs = []

    for env, heat, dhw, pv in itertools.product(
        env_options,
        heat_options,
        dhw_options,
        pv_options
    ):
        configs.append({
            "env": env,
            "heat": heat,
            "dhw": dhw,
            "fuel": b0["fuel"],
            "pv_percentage": pv
        })

    return configs


def optimize_configuration_per_building(
    config_current,
    current_dictionary,
    baseline_gdf_path,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20
):
    candidate_lists = {}
    indices = {}
    best_configs = {}
    best_npvs = {}

    for bid, b0 in config_current.items():
        cands = [b0] + list(generate_configs_for_building(b0))
        candidate_lists[bid] = cands
        indices[bid] = 0
        best_configs[bid] = b0
        best_npvs[bid] = 0.0

    active = True

    while active:
        active = False
        current_config = {}

        for bid in config_current:
            idx = indices[bid]
            cands = candidate_lists[bid]

            if idx < len(cands):
                current_config[bid] = cands[idx]
                active = True
            else:
                current_config[bid] = best_configs[bid]

        if not active:
            break

        result = compare_config_with_reference(
            reference_dictionary=current_dictionary,
            baseline_gdf_path=baseline_gdf_path,
            configuration=current_config,
            ee_measure_path=ee_measure_path,
            pv_type_path=pv_type_path,
            hp_catalog_path=hp_catalog_path,
            grid_pricing_path=grid_pricing_path,
            dhn_pricing_path=dhn_pricing_path,
            fuels_path=fuels_path,
            spot_price_path=spot_price_path,
            weatherfile_path=weatherfile_path,
            mycity=mycity,
            baseline_scenario=baseline_scenario,
            r=r,
            T=T
        )

        for bid in config_current:
            idx = indices[bid]
            cands = candidate_lists[bid]

            if idx < len(cands):
                conf = cands[idx]
                npv = result[bid]["costs"]["financial"]["NPV"]

                if npv > best_npvs[bid]:
                    best_npvs[bid] = npv
                    best_configs[bid] = conf

        for bid in config_current:
            if indices[bid] < len(candidate_lists[bid]):
                indices[bid] += 1

    final_result = compare_config_with_reference(
        reference_dictionary=current_dictionary,
        baseline_gdf_path=baseline_gdf_path,
        configuration=best_configs,
        ee_measure_path=ee_measure_path,
        pv_type_path=pv_type_path,
        hp_catalog_path=hp_catalog_path,
        grid_pricing_path=grid_pricing_path,
        dhn_pricing_path=dhn_pricing_path,
        fuels_path=fuels_path,
        spot_price_path=spot_price_path,
        weatherfile_path=weatherfile_path,
        mycity=mycity,
        baseline_scenario=baseline_scenario,
        r=r,
        T=T
    )

    return best_configs, final_result


def optimize_configuration_per_building_one_step(
    config_current,
    current_dictionary,
    baseline_dictionary,
    baseline_gdf_path,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20
):
    candidate_lists = {}
    indices = {}
    best_configs = {}
    best_npvs = {}

    for bid, b0 in config_current.items():
        cands = [b0] + list(generate_configs_for_building(b0))
        candidate_lists[bid] = cands
        indices[bid] = 0
        best_configs[bid] = b0
        best_npvs[bid] = 0.0

    active = True

    while active:
        active = False
        current_config = {}

        for bid in config_current:
            idx = indices[bid]
            cands = candidate_lists[bid]

            if idx < len(cands):
                current_config[bid] = cands[idx]
                active = True
            else:
                current_config[bid] = best_configs[bid]

        if not active:
            break

        result = compare_config_with_reference(
            reference_dictionary=current_dictionary,
            baseline_gdf_path=baseline_gdf_path,
            configuration=current_config,
            ee_measure_path=ee_measure_path,
            pv_type_path=pv_type_path,
            hp_catalog_path=hp_catalog_path,
            grid_pricing_path=grid_pricing_path,
            dhn_pricing_path=dhn_pricing_path,
            fuels_path=fuels_path,
            spot_price_path=spot_price_path,
            weatherfile_path=weatherfile_path,
            mycity=mycity,
            baseline_scenario=baseline_scenario,
            r=r,
            T=T
        )

        for bid in config_current:
            idx = indices[bid]
            cands = candidate_lists[bid]

            if idx < len(cands):
                conf = cands[idx]
                npv = result[bid]["costs"]["financial"]["NPV"]

                if npv > best_npvs[bid]:
                    best_npvs[bid] = npv
                    best_configs[bid] = conf

        for bid in config_current:
            if indices[bid] < len(candidate_lists[bid]):
                indices[bid] += 1

    baseline_gdf = gpd.read_file(baseline_gdf_path)

    id_col = "id"
    occ_col = "Number of occupants"

    occupants_map = (
        baseline_gdf[[id_col, occ_col]]
        .drop_duplicates(subset=[id_col])
        .set_index(id_col)[occ_col]
        .to_dict()
    )

    best_bid_global = None
    best_score_global = -math.inf

    for bid in config_current:
        n_occ = occupants_map.get(bid, None)

        if n_occ is None or n_occ <= 0:
            continue

        score = best_npvs[bid] / n_occ

        if score > best_score_global:
            best_score_global = score
            best_bid_global = bid

    final_configs = {}

    for bid, b0 in config_current.items():
        if bid == best_bid_global:
            final_configs[bid] = best_configs[bid]
        else:
            final_configs[bid] = b0

    final_result = compare_config_with_reference(
        reference_dictionary=baseline_dictionary,
        baseline_gdf_path=baseline_gdf_path,
        configuration=final_configs,
        ee_measure_path=ee_measure_path,
        pv_type_path=pv_type_path,
        hp_catalog_path=hp_catalog_path,
        grid_pricing_path=grid_pricing_path,
        dhn_pricing_path=dhn_pricing_path,
        fuels_path=fuels_path,
        spot_price_path=spot_price_path,
        weatherfile_path=weatherfile_path,
        mycity=mycity,
        baseline_scenario=baseline_scenario,
        r=r,
        T=T
    )

    return final_configs, final_result

import math
import copy
import pandas as pd
import geopandas as gpd


ENV_ORDER = ["none", "shallow", "medium", "deep"]
HP_ORDER = ["hp_le", "hp_me", "hp_he"]


def _conf_get(conf, names, default=None):
    for name in names:
        if name in conf:
            return conf[name]
    return default


def _conf_signature(conf):
    return (
        _conf_get(conf, ["env", "envelope", "envelope_measure"]),
        _conf_get(conf, ["sh_source", "SHsource", "SH_source", "space_heating_source"]),
        _conf_get(conf, ["dhw_source", "DHWsource", "DHW_source"]),
        _conf_get(conf, ["pv_percentage", "pv_percent", "PV_percentage", "pv"], 0),
    )


def _to_float(x, default=0.0):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def _env_rank(x):
    if x in ENV_ORDER:
        return ENV_ORDER.index(x)
    return 0


def _hp_rank(x):
    if x in HP_ORDER:
        return HP_ORDER.index(x) + 1
    return 0


def _is_hp(x):
    return isinstance(x, str) and x.startswith("hp")


def _proxy_score_candidate(b0, cand):
    env0, sh0, dhw0, pv0 = _conf_signature(b0)
    env1, sh1, dhw1, pv1 = _conf_signature(cand)

    score = 0.0

    env_gain = _env_rank(env1) - _env_rank(env0)
    if env_gain > 0:
        score += 1000.0 * env_gain
    elif env_gain < 0:
        score -= 100000.0

    pv_gain = _to_float(pv1) - _to_float(pv0)
    if pv_gain > 0:
        score += 15.0 * pv_gain
    elif pv_gain < 0:
        score -= 100000.0

    for old, new in [(sh0, sh1), (dhw0, dhw1)]:
        if old == new:
            score += 0.0
        elif new == "dhn":
            score += 700.0
        elif _is_hp(new):
            score += 600.0 + 100.0 * _hp_rank(new)
        elif new == "boiler":
            score -= 200.0

        if old == "boiler" and new != "boiler":
            score += 800.0

        if _is_hp(old) and _is_hp(new):
            score += 100.0 * (_hp_rank(new) - _hp_rank(old))

    changed = 0
    changed += env0 != env1
    changed += sh0 != sh1
    changed += dhw0 != dhw1
    changed += _to_float(pv0) != _to_float(pv1)

    score -= 25.0 * changed

    return score


def reduce_candidates_for_building(b0, candidates, top_k=5, keep_diverse=True):
    current_sig = _conf_signature(b0)

    unique = {}
    for cand in candidates:
        sig = _conf_signature(cand)
        if sig not in unique:
            unique[sig] = cand

    scored = []

    for sig, cand in unique.items():
        if sig == current_sig:
            continue

        scored.append(
            {
                "candidate": cand,
                "signature": sig,
                "proxy_score": _proxy_score_candidate(b0, cand),
            }
        )

    scored = sorted(scored, key=lambda x: x["proxy_score"], reverse=True)

    if not keep_diverse:
        return [b0] + [x["candidate"] for x in scored[:top_k]]

    selected = []
    used_env = set()
    used_heat = set()
    used_pv = set()

    for item in scored:
        env, sh, dhw, pv = item["signature"]
        heat_sig = (sh, dhw)

        diverse = (
            env not in used_env
            or heat_sig not in used_heat
            or pv not in used_pv
        )

        if diverse or len(selected) < max(1, top_k // 2):
            selected.append(item["candidate"])
            used_env.add(env)
            used_heat.add(heat_sig)
            used_pv.add(pv)

        if len(selected) >= top_k:
            break

    if len(selected) < top_k:
        for item in scored:
            cand = item["candidate"]
            if cand not in selected:
                selected.append(cand)
            if len(selected) >= top_k:
                break

    return [b0] + selected


def load_occupants_and_area_maps(baseline_gdf_path, score_mode, id_col="id"):
    baseline_gdf = gpd.read_file(baseline_gdf_path)

    occupants_map = {}

    if "Number of occupants" in baseline_gdf.columns:
        occupants_map = (
            baseline_gdf[[id_col, "Number of occupants"]]
            .drop_duplicates(subset=[id_col])
            .set_index(id_col)["Number of occupants"]
            .to_dict()
        )

    area_map = None

    if score_mode == "npv_per_m2":
        area_gdf = baseline_gdf.copy()

        if area_gdf.crs is None:
            area_gdf = area_gdf.set_crs("EPSG:3006")

        if area_gdf.crs.is_geographic:
            area_calc = area_gdf.to_crs("EPSG:3006")
        else:
            area_calc = area_gdf

        area_gdf["footprint_area_m2"] = area_calc.geometry.area

        if "Floors" in area_gdf.columns:
            floors = pd.to_numeric(area_gdf["Floors"], errors="coerce").fillna(1.0)
        else:
            floors = pd.Series(1.0, index=area_gdf.index)

        floors = floors.clip(lower=1.0)
        area_gdf["total_floor_area_m2"] = area_gdf["footprint_area_m2"] * floors

        area_map = (
            area_gdf[[id_col, "total_floor_area_m2"]]
            .drop_duplicates(subset=[id_col])
            .set_index(id_col)["total_floor_area_m2"]
            .to_dict()
        )

    return occupants_map, area_map


def get_financial_npv(result, bid, default=0.0):
    try:
        return result[bid]["costs"]["financial"]["NPV"]
    except Exception:
        return default


def score_mover(bid, improvement_npv, score_mode, occupants_map=None, area_map=None):
    if score_mode == "npv":
        return improvement_npv

    if score_mode == "npv_per_occupant":
        n_occ = occupants_map.get(bid, None) if occupants_map is not None else None
        n_occ = _to_float(n_occ, default=0.0)

        if n_occ <= 0:
            return None

        return improvement_npv / n_occ

    if score_mode == "npv_per_m2":
        area = area_map.get(bid, None) if area_map is not None else None
        area = _to_float(area, default=0.0)

        if area <= 0:
            return None

        return improvement_npv / area

    raise ValueError(
        "score_mode must be one of: 'npv', 'npv_per_occupant', 'npv_per_m2'"
    )


def select_movers(
    config_current,
    best_configs,
    best_npvs,
    baseline_gdf_path,
    move_fraction=0.10,
    min_movers=1,
    score_mode="npv_per_occupant",
):
    occupants_map, area_map = load_occupants_and_area_maps(
        baseline_gdf_path=baseline_gdf_path,
        score_mode=score_mode,
    )

    movers = []

    for bid in config_current:
        improvement_npv = best_npvs.get(bid, 0.0)

        if improvement_npv <= 0:
            continue

        if best_configs.get(bid) == config_current[bid]:
            continue

        score = score_mover(
            bid=bid,
            improvement_npv=improvement_npv,
            score_mode=score_mode,
            occupants_map=occupants_map,
            area_map=area_map,
        )

        if score is None:
            continue

        movers.append(
            {
                "bid": bid,
                "best_config": best_configs[bid],
                "improvement_npv": improvement_npv,
                "score": score,
            }
        )

    movers = sorted(movers, key=lambda x: x["score"], reverse=True)

    if not movers:
        return [], set()

    n_to_move = math.ceil(len(movers) * move_fraction)
    n_to_move = max(min_movers, n_to_move)
    n_to_move = min(n_to_move, len(movers))

    selected_movers = movers[:n_to_move]
    selected_bids = {m["bid"] for m in selected_movers}

    return selected_movers, selected_bids


def add_step_metadata(final_result, best_npvs, selected_movers, selected_bids):
    selected_lookup = {m["bid"]: m for m in selected_movers}

    for bid, data in final_result.items():
        if "costs" not in data:
            data["costs"] = {}

        if "financial" not in data["costs"]:
            data["costs"]["financial"] = {}

        data["costs"]["financial"]["step_best_improvement_npv"] = best_npvs.get(bid, 0.0)
        data["costs"]["financial"]["step_selected_to_move"] = bid in selected_bids

        if bid in selected_lookup:
            data["costs"]["financial"]["step_selection_score"] = selected_lookup[bid]["score"]
        else:
            data["costs"]["financial"]["step_selection_score"] = 0.0

    return final_result

import math
import copy
import pandas as pd
import geopandas as gpd


ENV_ORDER = ["none", "shallow", "medium", "deep"]

ENV_ALLOWED = {
    "none": ["none", "shallow", "medium", "deep"],
    "shallow": ["shallow", "medium", "deep"],
    "medium": ["medium", "deep"],
    "deep": ["deep"],
}

HEAT_ALLOWED = {
    "boiler": ["boiler", "dhn", "hp_le", "hp_me", "hp_he"],
    "dhn": ["dhn", "hp_le", "hp_me", "hp_he"],
    "hp_le": ["hp_le", "hp_me", "hp_he", "dhn"],
    "hp_me": ["hp_me", "hp_he", "dhn"],
    "hp_he": ["hp_he", "dhn"],
    None: ["dhn", "hp_le", "hp_me", "hp_he", "boiler"],
}

SEARCH_STAGES = [
    "env",
    "heat",
    "pv_coarse",
    "pv_fine",
]

def evaluate_candidate_pool_by_step(
    config_current,
    candidate_pool_by_building,
    current_dictionary,
    baseline_gdf_path,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20,
    verbose=True,
):
    max_len = max(len(v) for v in candidate_pool_by_building.values())

    evaluated = {bid: [] for bid in candidate_pool_by_building}

    for step in range(max_len):
        test_config = {}

        for bid, b0 in config_current.items():
            pool = candidate_pool_by_building.get(bid, [b0])

            if step < len(pool):
                test_config[bid] = pool[step]
            else:
                test_config[bid] = pool[-1]

        if verbose:
            print("evaluation step", step + 1, "of", max_len)

        result = compare_config_with_reference(
            reference_dictionary=current_dictionary,
            baseline_gdf_path=baseline_gdf_path,
            configuration=test_config,
            ee_measure_path=ee_measure_path,
            pv_type_path=pv_type_path,
            hp_catalog_path=hp_catalog_path,
            grid_pricing_path=grid_pricing_path,
            dhn_pricing_path=dhn_pricing_path,
            fuels_path=fuels_path,
            spot_price_path=spot_price_path,
            weatherfile_path=weatherfile_path,
            mycity=mycity,
            baseline_scenario=baseline_scenario,
            r=r,
            T=T,
        )

        for bid in candidate_pool_by_building:
            pool = candidate_pool_by_building[bid]

            if step >= len(pool):
                continue

            conf = pool[step]

            try:
                npv = result[bid]["costs"]["financial"]["NPV"]
            except Exception:
                npv = -math.inf

            evaluated[bid].append(
                {
                    "config": conf,
                    "npv": npv,
                }
            )

    return evaluated

def run_beam_search_candidates(
    config_current,
    current_dictionary,
    baseline_gdf_path,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20,
    beam_width=5,
    search_stages=None,
    verbose=True,
):
    if search_stages is None:
        search_stages = ["env", "heat", "pv_coarse", "pv_fine"]

    beams = {
        bid: [
            {
                "config": b0,
                "npv": 0.0,
            }
        ]
        for bid, b0 in config_current.items()
    }

    for stage in search_stages:
        if verbose:
            print("beam stage:", stage)

        candidate_pool_by_building = {}

        for bid, items in beams.items():
            expanded = []

            for item in items:
                expanded.extend(generate_stage_candidates(item["config"], stage))

            expanded = unique_configs(expanded)

            candidate_pool_by_building[bid] = expanded

        evaluated = evaluate_candidate_pool_by_step(
            config_current=config_current,
            candidate_pool_by_building=candidate_pool_by_building,
            current_dictionary=current_dictionary,
            baseline_gdf_path=baseline_gdf_path,
            ee_measure_path=ee_measure_path,
            pv_type_path=pv_type_path,
            hp_catalog_path=hp_catalog_path,
            grid_pricing_path=grid_pricing_path,
            dhn_pricing_path=dhn_pricing_path,
            fuels_path=fuels_path,
            spot_price_path=spot_price_path,
            weatherfile_path=weatherfile_path,
            mycity=mycity,
            baseline_scenario=baseline_scenario,
            r=r,
            T=T,
            verbose=verbose,
        )

        new_beams = {}

        for bid, items in evaluated.items():
            items = sorted(items, key=lambda x: x["npv"], reverse=True)
            new_beams[bid] = items[:beam_width]

        beams = new_beams

    best_configs = {}
    best_npvs = {}

    for bid, items in beams.items():
        items = sorted(items, key=lambda x: x["npv"], reverse=True)
        best = items[0]

        if best["npv"] > 0:
            best_configs[bid] = best["config"]
            best_npvs[bid] = best["npv"]
        else:
            best_configs[bid] = config_current[bid]
            best_npvs[bid] = 0.0

    return best_configs, best_npvs, beams

def load_occupants_and_area_maps(baseline_gdf_path, score_mode, id_col="id"):
    baseline_gdf = gpd.read_file(baseline_gdf_path)

    occupants_map = {}

    if "Number of occupants" in baseline_gdf.columns:
        occupants_map = (
            baseline_gdf[[id_col, "Number of occupants"]]
            .drop_duplicates(subset=[id_col])
            .set_index(id_col)["Number of occupants"]
            .to_dict()
        )

    area_map = None

    if score_mode == "npv_per_m2":
        area_gdf = baseline_gdf.copy()

        if area_gdf.crs is None:
            area_gdf = area_gdf.set_crs("EPSG:3006")

        if area_gdf.crs.is_geographic:
            area_calc = area_gdf.to_crs("EPSG:3006")
        else:
            area_calc = area_gdf

        area_gdf["footprint_area_m2"] = area_calc.geometry.area

        if "Floors" in area_gdf.columns:
            floors = pd.to_numeric(area_gdf["Floors"], errors="coerce").fillna(1.0)
        else:
            floors = pd.Series(1.0, index=area_gdf.index)

        floors = floors.clip(lower=1.0)
        area_gdf["total_floor_area_m2"] = area_gdf["footprint_area_m2"] * floors

        area_map = (
            area_gdf[[id_col, "total_floor_area_m2"]]
            .drop_duplicates(subset=[id_col])
            .set_index(id_col)["total_floor_area_m2"]
            .to_dict()
        )

    return occupants_map, area_map


def safe_float(x, default=0.0):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def calculate_selection_score(
    bid,
    improvement_npv,
    score_mode,
    occupants_map=None,
    area_map=None,
):
    if score_mode == "npv":
        return improvement_npv

    if score_mode == "npv_per_occupant":
        n_occ = safe_float(occupants_map.get(bid, 0.0), 0.0)

        if n_occ <= 0:
            return None

        return improvement_npv / n_occ

    if score_mode == "npv_per_m2":
        area = safe_float(area_map.get(bid, 0.0), 0.0)

        if area <= 0:
            return None

        return improvement_npv / area

    raise ValueError(
        "score_mode must be one of: 'npv', 'npv_per_occupant', 'npv_per_m2'"
    )


def select_movers_from_best_configs(
    config_current,
    best_configs,
    best_npvs,
    baseline_gdf_path,
    move_fraction=0.10,
    min_movers=1,
    score_mode="npv_per_occupant",
):
    occupants_map, area_map = load_occupants_and_area_maps(
        baseline_gdf_path=baseline_gdf_path,
        score_mode=score_mode,
    )

    movers = []

    for bid in config_current:
        improvement_npv = best_npvs.get(bid, 0.0)

        if improvement_npv <= 0:
            continue

        if best_configs[bid] == config_current[bid]:
            continue

        score = calculate_selection_score(
            bid=bid,
            improvement_npv=improvement_npv,
            score_mode=score_mode,
            occupants_map=occupants_map,
            area_map=area_map,
        )

        if score is None:
            continue

        movers.append(
            {
                "bid": bid,
                "best_config": best_configs[bid],
                "improvement_npv": improvement_npv,
                "score": score,
            }
        )

    movers = sorted(movers, key=lambda x: x["score"], reverse=True)

    if not movers:
        return [], set()

    n_to_move = math.ceil(len(movers) * move_fraction)
    n_to_move = max(min_movers, n_to_move)
    n_to_move = min(n_to_move, len(movers))

    selected_movers = movers[:n_to_move]
    selected_bids = {m["bid"] for m in selected_movers}

    return selected_movers, selected_bids


def add_step_metadata(final_result, best_npvs, selected_movers, selected_bids):
    selected_lookup = {m["bid"]: m for m in selected_movers}

    for bid, data in final_result.items():
        if "costs" not in data:
            data["costs"] = {}

        if "financial" not in data["costs"]:
            data["costs"]["financial"] = {}

        data["costs"]["financial"]["step_best_improvement_npv"] = best_npvs.get(bid, 0.0)
        data["costs"]["financial"]["step_selected_to_move"] = bid in selected_bids

        if bid in selected_lookup:
            data["costs"]["financial"]["step_selection_score"] = selected_lookup[bid]["score"]
        else:
            data["costs"]["financial"]["step_selection_score"] = 0.0

    return final_result

def get_conf_value(conf, possible_keys, default=None):
    for key in possible_keys:
        if key in conf:
            return conf[key]
    return default


def set_conf_value(conf, possible_keys, value):
    out = copy.deepcopy(conf)

    for key in possible_keys:
        if key in out:
            out[key] = value
            return out

    out[possible_keys[0]] = value
    return out


def config_signature(conf):
    return (
        get_conf_value(conf, ["env", "envelope", "envelope_measure"]),
        get_conf_value(conf, ["sh_source", "SHsource", "SH_source", "space_heating_source"]),
        get_conf_value(conf, ["dhw_source", "DHWsource", "DHW_source"]),
        get_conf_value(conf, ["pv_percentage", "pv_percent", "PV_percentage"], 0),
    )


def unique_configs(configs):
    seen = set()
    out = []

    for conf in configs:
        sig = config_signature(conf)

        if sig not in seen:
            seen.add(sig)
            out.append(conf)

    return out


def generate_stage_candidates(conf, stage, pv_fine_step=10):
    env, sh, dhw, pv = config_signature(conf)

    candidates = []

    if stage == "env":
        allowed_envs = ENV_ALLOWED.get(env, ENV_ALLOWED["none"])

        for new_env in allowed_envs:
            c = set_conf_value(conf, ["env", "envelope", "envelope_measure"], new_env)
            candidates.append(c)

    elif stage == "heat":
        allowed_sh = HEAT_ALLOWED.get(sh, HEAT_ALLOWED[None])
        allowed_dhw = HEAT_ALLOWED.get(dhw, HEAT_ALLOWED[None])

        for new_sh in allowed_sh:
            for new_dhw in allowed_dhw:
                c = set_conf_value(conf, ["sh_source", "SHsource", "SH_source", "space_heating_source"], new_sh)
                c = set_conf_value(c, ["dhw_source", "DHWsource", "DHW_source"], new_dhw)
                candidates.append(c)

    elif stage == "pv_coarse":
        current_pv = float(pv or 0)

        for new_pv in [current_pv, 0, 25, 50, 75, 100]:
            if new_pv >= current_pv:
                c = set_conf_value(conf, ["pv_percentage", "pv_percent", "PV_percentage"], new_pv)
                candidates.append(c)

    elif stage == "pv_fine":
        current_pv = float(pv or 0)

        values = [
            current_pv - 2 * pv_fine_step,
            current_pv - pv_fine_step,
            current_pv,
            current_pv + pv_fine_step,
            current_pv + 2 * pv_fine_step,
        ]

        values = [v for v in values if 0 <= v <= 100]

        for new_pv in values:
            c = set_conf_value(conf, ["pv_percentage", "pv_percent", "PV_percentage"], new_pv)
            candidates.append(c)

    else:
        raise ValueError(f"Unknown search stage: {stage}")

    return unique_configs(candidates)


def optimize_configuration_per_building_batch_step(
    config_current,
    current_dictionary,
    baseline_dictionary,
    baseline_gdf_path,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20,
    move_fraction=0.10,
    min_movers=1,
    score_mode="npv_per_occupant",
    pv_step=10,
    allow_deep_env_jump=False,
    max_local_iterations=20,
    min_improvement=0.0,
    verbose=True,
):
    if move_fraction <= 0:
        raise ValueError("move_fraction must be > 0")

    if move_fraction > 1:
        raise ValueError("move_fraction must be <= 1")

    converged_configs, best_configs, best_npvs, history = local_best_response_until_convergence(
        config_current=config_current,
        current_dictionary=current_dictionary,
        baseline_gdf_path=baseline_gdf_path,
        ee_measure_path=ee_measure_path,
        pv_type_path=pv_type_path,
        hp_catalog_path=hp_catalog_path,
        grid_pricing_path=grid_pricing_path,
        dhn_pricing_path=dhn_pricing_path,
        fuels_path=fuels_path,
        spot_price_path=spot_price_path,
        weatherfile_path=weatherfile_path,
        mycity=mycity,
        baseline_scenario=baseline_scenario,
        r=r,
        T=T,
        pv_step=pv_step,
        allow_deep_env_jump=allow_deep_env_jump,
        max_local_iterations=max_local_iterations,
        min_improvement=min_improvement,
        verbose=verbose,
    )

    converged_result = compare_config_with_reference(
        reference_dictionary=current_dictionary,
        baseline_gdf_path=baseline_gdf_path,
        configuration=converged_configs,
        ee_measure_path=ee_measure_path,
        pv_type_path=pv_type_path,
        hp_catalog_path=hp_catalog_path,
        grid_pricing_path=grid_pricing_path,
        dhn_pricing_path=dhn_pricing_path,
        fuels_path=fuels_path,
        spot_price_path=spot_price_path,
        weatherfile_path=weatherfile_path,
        mycity=mycity,
        baseline_scenario=baseline_scenario,
        r=r,
        T=T,
    )

    selected_movers, selected_bids = select_top_movers_after_convergence(
        config_current=config_current,
        converged_configs=converged_configs,
        converged_result=converged_result,
        baseline_gdf_path=baseline_gdf_path,
        move_fraction=move_fraction,
        min_movers=min_movers,
        score_mode=score_mode,
    )

    final_configs = {}

    for bid, b0 in config_current.items():
        if bid in selected_bids:
            final_configs[bid] = converged_configs[bid]
        else:
            final_configs[bid] = b0

    if verbose:
        print("local-search iterations:", len(history))
        print("candidate movers:", len(selected_movers))
        print("selected movers:", len(selected_bids))

    final_result = compare_config_with_reference(
        reference_dictionary=baseline_dictionary,
        baseline_gdf_path=baseline_gdf_path,
        configuration=final_configs,
        ee_measure_path=ee_measure_path,
        pv_type_path=pv_type_path,
        hp_catalog_path=hp_catalog_path,
        grid_pricing_path=grid_pricing_path,
        dhn_pricing_path=dhn_pricing_path,
        fuels_path=fuels_path,
        spot_price_path=spot_price_path,
        weatherfile_path=weatherfile_path,
        mycity=mycity,
        baseline_scenario=baseline_scenario,
        r=r,
        T=T,
    )

    final_result = add_step_metadata(
        final_result=final_result,
        converged_result=converged_result,
        selected_movers=selected_movers,
        selected_bids=selected_bids,
    )

    return final_configs, final_result, history
import copy
import math
import pandas as pd
import geopandas as gpd


ENV_ORDER = ["none", "shallow", "medium", "deep"]
HP_ORDER = ["hp_le", "hp_me", "hp_he"]


def get_conf_value(conf, keys, default=None):
    for key in keys:
        if key in conf:
            return conf[key]
    return default


def set_conf_value(conf, keys, value):
    out = copy.deepcopy(conf)

    for key in keys:
        if key in out:
            out[key] = value
            return out

    out[keys[0]] = value
    return out


def config_signature(conf):
    return (
        get_conf_value(conf, ["env", "envelope", "envelope_measure"]),
        get_conf_value(conf, ["sh_source", "SHsource", "SH_source", "space_heating_source"]),
        get_conf_value(conf, ["dhw_source", "DHWsource", "DHW_source"]),
        float(get_conf_value(conf, ["pv_percentage", "pv_percent", "PV_percentage"], 0) or 0),
    )


def configs_equal(a, b):
    return config_signature(a) == config_signature(b)


def unique_configs(configs):
    seen = set()
    out = []

    for conf in configs:
        sig = config_signature(conf)

        if sig not in seen:
            seen.add(sig)
            out.append(conf)

    return out


def next_env_values(env, allow_deep_env_jump=False):
    if env not in ENV_ORDER:
        env = "none"

    i = ENV_ORDER.index(env)
    values = []

    if i + 1 < len(ENV_ORDER):
        values.append(ENV_ORDER[i + 1])

    if allow_deep_env_jump:
        values.extend(ENV_ORDER[i + 2:])

    return values


def heat_neighbor_values(source):
    values = []

    if source == "boiler":
        values.extend(["dhn", "hp_le"])

    elif source == "dhn":
        values.append("hp_le")

    elif source in HP_ORDER:
        i = HP_ORDER.index(source)

        if i + 1 < len(HP_ORDER):
            values.append(HP_ORDER[i + 1])

        values.append("dhn")

    else:
        values.extend(["dhn", "hp_le", "boiler"])

    return list(dict.fromkeys(values))




ENV_KEY = "env"
HEAT_KEY = "heat"
DHW_KEY = "dhw"
PV_KEY = "pv_percentage"

ENV_ORDER = ["none", "shallow", "medium", "deep"]
HEAT_TARGETS = ["dhn", "hp_le", "hp_me", "hp_he"]


def config_signature(conf):
    return (
        conf[ENV_KEY],
        conf[HEAT_KEY],
        conf[DHW_KEY],
        float(conf.get(PV_KEY, 0.0) or 0.0),
    )


def configs_equal(a, b):
    return config_signature(a) == config_signature(b)


def set_config_value(conf, key, value):
    if key not in conf:
        raise KeyError(f"Missing config key: {key}. Available keys: {list(conf.keys())}")

    out = copy.deepcopy(conf)
    out[key] = value
    return out


def unique_configs(configs):
    seen = set()
    out = []

    for conf in configs:
        sig = config_signature(conf)

        if sig not in seen:
            seen.add(sig)
            out.append(conf)

    return out


def generate_neighbor_configs(
    conf,
    pv_step=10,
    include_current=True,
    include_packages=True,
):
    env, heat, dhw, pv = config_signature(conf)

    if env not in ENV_ORDER:
        env = "none"

    env_i = ENV_ORDER.index(env)

    neighbors = []

    if include_current:
        neighbors.append(copy.deepcopy(conf))

    # Envelope jumps: allow all upgrades, not only +1.
    for new_env in ENV_ORDER[env_i + 1:]:
        c = set_config_value(conf, ENV_KEY, new_env)
        neighbors.append(c)

    # Heat-only changes.
    for target in HEAT_TARGETS:
        if target != heat:
            c = set_config_value(conf, HEAT_KEY, target)
            neighbors.append(c)

    # DHW-only changes.
    for target in HEAT_TARGETS:
        if target != dhw:
            c = set_config_value(conf, DHW_KEY, target)
            neighbors.append(c)

    # Heat + DHW package changes.
    for target in HEAT_TARGETS:
        if target != heat or target != dhw:
            c = set_config_value(conf, HEAT_KEY, target)
            c = set_config_value(c, DHW_KEY, target)
            neighbors.append(c)

    # PV jumps.
    pv_levels = [
        pv + pv_step,
        25.0,
        50.0,
        75.0,
        100.0,
    ]

    for new_pv in pv_levels:
        new_pv = float(new_pv)

        if pv <= new_pv <= 100.0:
            c = set_config_value(conf, PV_KEY, new_pv)
            neighbors.append(c)

    if include_packages:
        packages = [
            ("shallow", "dhn", "dhn", 25.0),
            ("medium", "dhn", "dhn", 50.0),
            ("deep", "dhn", "dhn", 75.0),
            ("deep", "dhn", "dhn", 100.0),

            ("shallow", "hp_le", "hp_le", 25.0),
            ("medium", "hp_me", "hp_me", 50.0),
            ("deep", "hp_he", "hp_he", 75.0),
            ("deep", "hp_he", "hp_he", 100.0),

            ("medium", "hp_me", "dhn", 50.0),
            ("medium", "dhn", "hp_me", 50.0),
            ("deep", "hp_he", "dhn", 75.0),
            ("deep", "dhn", "hp_he", 75.0),
        ]

        for new_env, new_heat, new_dhw, new_pv in packages:
            if ENV_ORDER.index(new_env) < env_i:
                continue

            if float(new_pv) < pv:
                continue

            c = set_config_value(conf, ENV_KEY, new_env)
            c = set_config_value(c, HEAT_KEY, new_heat)
            c = set_config_value(c, DHW_KEY, new_dhw)
            c = set_config_value(c, PV_KEY, float(new_pv))
            neighbors.append(c)

    return unique_configs(neighbors)


def evaluate_neighbor_pools(
    working_configs,
    neighbor_pools,
    current_dictionary,
    baseline_gdf_path,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20,
    verbose=True,
):
    max_neighbors = max(len(pool) for pool in neighbor_pools.values())

    evaluated = {bid: [] for bid in working_configs}

    for step in range(max_neighbors):
        test_config = {}

        for bid, current_conf in working_configs.items():
            pool = neighbor_pools[bid]

            if step < len(pool):
                test_config[bid] = pool[step]
            else:
                test_config[bid] = current_conf

        if verbose:
            print("neighbor evaluation", step + 1, "of", max_neighbors)

        result = compare_config_with_reference(
            reference_dictionary=current_dictionary,
            baseline_gdf_path=baseline_gdf_path,
            configuration=test_config,
            ee_measure_path=ee_measure_path,
            pv_type_path=pv_type_path,
            hp_catalog_path=hp_catalog_path,
            grid_pricing_path=grid_pricing_path,
            dhn_pricing_path=dhn_pricing_path,
            fuels_path=fuels_path,
            spot_price_path=spot_price_path,
            weatherfile_path=weatherfile_path,
            mycity=mycity,
            baseline_scenario=baseline_scenario,
            r=r,
            T=T,
        )

        for bid in working_configs:
            pool = neighbor_pools[bid]

            if step >= len(pool):
                continue

            try:
                npv = result[bid]["costs"]["financial"]["NPV"]
            except Exception:
                npv = -math.inf

            evaluated[bid].append(
                {
                    "config": pool[step],
                    "npv": npv,
                }
            )

    return evaluated


def local_best_response_until_convergence(
    config_current,
    current_dictionary,
    baseline_gdf_path,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path,
    weatherfile_path,
    mycity,
    baseline_scenario,
    r=0.04,
    T=20,
    pv_step=10,
    allow_deep_env_jump=False,
    max_local_iterations=20,
    min_improvement=0.0,
    verbose=True,
):
    working_configs = copy.deepcopy(config_current)

    best_configs = copy.deepcopy(config_current)
    best_npvs = {bid: 0.0 for bid in config_current}

    history = []

    for iteration in range(max_local_iterations):
        if verbose:
            print("local best-response iteration", iteration + 1, "of", max_local_iterations)

        neighbor_pools = {}

        for bid, conf in working_configs.items():
            neighbor_pools[bid] = generate_neighbor_configs(
                conf=conf,
                pv_step=pv_step,
                include_current=True,
            )

        evaluated = evaluate_neighbor_pools(
            working_configs=working_configs,
            neighbor_pools=neighbor_pools,
            current_dictionary=current_dictionary,
            baseline_gdf_path=baseline_gdf_path,
            ee_measure_path=ee_measure_path,
            pv_type_path=pv_type_path,
            hp_catalog_path=hp_catalog_path,
            grid_pricing_path=grid_pricing_path,
            dhn_pricing_path=dhn_pricing_path,
            fuels_path=fuels_path,
            spot_price_path=spot_price_path,
            weatherfile_path=weatherfile_path,
            mycity=mycity,
            baseline_scenario=baseline_scenario,
            r=r,
            T=T,
            verbose=verbose,
        )

        next_configs = copy.deepcopy(working_configs)
        changed_bids = []

        for bid, items in evaluated.items():
            if not items:
                continue

            current_item = None

            for item in items:
                if configs_equal(item["config"], working_configs[bid]):
                    current_item = item
                    break

            if current_item is None:
                current_npv = 0.0
            else:
                current_npv = current_item["npv"]

            best_item = max(items, key=lambda x: x["npv"])

            if best_item["npv"] > current_npv + min_improvement:
                next_configs[bid] = best_item["config"]
                changed_bids.append(bid)

            if best_item["npv"] > best_npvs[bid]:
                best_npvs[bid] = best_item["npv"]
                best_configs[bid] = best_item["config"]

        history.append(
            {
                "iteration": iteration + 1,
                "changed_buildings": len(changed_bids),
                "changed_bids": changed_bids,
            }
        )

        if verbose:
            print("changed buildings:", len(changed_bids))

        if len(changed_bids) == 0:
            break

        working_configs = next_configs

    return working_configs, best_configs, best_npvs, history

def load_occupants_and_area_maps(baseline_gdf_path, score_mode, id_col="id"):
    baseline_gdf = gpd.read_file(baseline_gdf_path)

    occupants_map = {}

    if "Number of occupants" in baseline_gdf.columns:
        occupants_map = (
            baseline_gdf[[id_col, "Number of occupants"]]
            .drop_duplicates(subset=[id_col])
            .set_index(id_col)["Number of occupants"]
            .to_dict()
        )

    area_map = None

    if score_mode == "npv_per_m2":
        area_gdf = baseline_gdf.copy()

        if area_gdf.crs is None:
            area_gdf = area_gdf.set_crs("EPSG:3006")

        if area_gdf.crs.is_geographic:
            area_calc = area_gdf.to_crs("EPSG:3006")
        else:
            area_calc = area_gdf

        area_gdf["footprint_area_m2"] = area_calc.geometry.area

        if "Floors" in area_gdf.columns:
            floors = pd.to_numeric(area_gdf["Floors"], errors="coerce").fillna(1.0)
        else:
            floors = pd.Series(1.0, index=area_gdf.index)

        floors = floors.clip(lower=1.0)
        area_gdf["total_floor_area_m2"] = area_gdf["footprint_area_m2"] * floors

        area_map = (
            area_gdf[[id_col, "total_floor_area_m2"]]
            .drop_duplicates(subset=[id_col])
            .set_index(id_col)["total_floor_area_m2"]
            .to_dict()
        )

    return occupants_map, area_map


def safe_float(x, default=0.0):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def calculate_selection_score(
    bid,
    improvement_npv,
    score_mode,
    occupants_map=None,
    area_map=None,
):
    if score_mode == "npv":
        return improvement_npv

    if score_mode == "npv_per_occupant":
        n_occ = safe_float(occupants_map.get(bid, 0.0), 0.0)

        if n_occ <= 0:
            return None

        return improvement_npv / n_occ

    if score_mode == "npv_per_m2":
        area = safe_float(area_map.get(bid, 0.0), 0.0)

        if area <= 0:
            return None

        return improvement_npv / area

    raise ValueError(
        "score_mode must be one of: 'npv', 'npv_per_occupant', 'npv_per_m2'"
    )


def select_top_movers_after_convergence(
    config_current,
    converged_configs,
    converged_result,
    baseline_gdf_path,
    move_fraction=0.10,
    min_movers=1,
    score_mode="npv_per_occupant",
):
    occupants_map, area_map = load_occupants_and_area_maps(
        baseline_gdf_path=baseline_gdf_path,
        score_mode=score_mode,
    )

    movers = []

    for bid in config_current:
        if configs_equal(config_current[bid], converged_configs[bid]):
            continue

        try:
            improvement_npv = converged_result[bid]["costs"]["financial"]["NPV"]
        except Exception:
            improvement_npv = 0.0

        if improvement_npv <= 0:
            continue

        score = calculate_selection_score(
            bid=bid,
            improvement_npv=improvement_npv,
            score_mode=score_mode,
            occupants_map=occupants_map,
            area_map=area_map,
        )

        if score is None:
            continue

        movers.append(
            {
                "bid": bid,
                "best_config": converged_configs[bid],
                "improvement_npv": improvement_npv,
                "score": score,
            }
        )

    movers = sorted(movers, key=lambda x: x["score"], reverse=True)

    if not movers:
        return [], set()

    n_to_move = math.ceil(len(movers) * move_fraction)
    n_to_move = max(min_movers, n_to_move)
    n_to_move = min(n_to_move, len(movers))

    selected_movers = movers[:n_to_move]
    selected_bids = {m["bid"] for m in selected_movers}

    return selected_movers, selected_bids


def add_step_metadata(final_result, converged_result, selected_movers, selected_bids):
    selected_lookup = {m["bid"]: m for m in selected_movers}

    for bid, data in final_result.items():
        if "costs" not in data:
            data["costs"] = {}

        if "financial" not in data["costs"]:
            data["costs"]["financial"] = {}

        try:
            convergence_npv = converged_result[bid]["costs"]["financial"]["NPV"]
        except Exception:
            convergence_npv = 0.0

        data["costs"]["financial"]["step_converged_improvement_npv"] = convergence_npv
        data["costs"]["financial"]["step_selected_to_move"] = bid in selected_bids

        if bid in selected_lookup:
            data["costs"]["financial"]["step_selection_score"] = selected_lookup[bid]["score"]
        else:
            data["costs"]["financial"]["step_selection_score"] = 0.0

    return final_result




def extract_config(gdf):
    def get_fuel(val):
        if isinstance(val, str) and "_" in val and "hp" not in val:
            return val.split("_", 1)[1]
        return None

    return {
        row["id"]: {
            "env": row["EEdepth"],
            "heat": row["SHSource"],
            "dhw": row["DHWsource"],
            "fuel": get_fuel(row["SHSource"]) or get_fuel(row["DHWsource"]),
            "pv_percentage": row["PVpercentage"]
        }
        for _, row in gdf.iterrows()
    }


def load_buildings(input_data):
    if isinstance(input_data, gpd.GeoDataFrame):
        gdf = input_data.copy()
    elif isinstance(input_data, str):
        gdf = gpd.read_file(input_data)
    else:
        raise ValueError("invalid_input")

    required_columns = [
        "Name",
        "id",
        "EEdepth",
        "SHSource",
        "DHWsource",
        "PVType",
        "PVpercentage"
    ]

    missing = [c for c in required_columns if c not in gdf.columns]

    if missing:
        raise ValueError(f"{missing}")

    return gdf


def p95(x):
    return np.percentile(x, 95)


def pick_hp(capacity, catalog):
    for hp in sorted(catalog, key=lambda x: x["max_heating_capacity_kW"]):
        if hp["max_heating_capacity_kW"] >= capacity:
            return hp

    return max(catalog, key=lambda x: x["max_heating_capacity_kW"])


def hp_cost(building, interv, hp_catalog, factor=1.0):
    base = building["base"]
    meta = building["meta"]

    sh_arr = np.array(base["space_heating"])
    dhw_arr = np.array(base["dhw_demand"])

    sh_before = interv.get("sh_source", (meta["SHSource"], None))[0]
    dhw_before = interv.get("dhw_source", (meta["DHWsource"], None))[0]

    sh_after = interv.get("sh_source", (None, None))[1]
    dhw_after = interv.get("dhw_source", (None, None))[1]

    sh_final = sh_after if sh_after is not None else meta["SHSource"]
    dhw_final = dhw_after if dhw_after is not None else meta["DHWsource"]

    sh_initial = sh_before if sh_before is not None else meta["SHSource"]
    dhw_initial = dhw_before if dhw_before is not None else meta["DHWsource"]

    sh_hp_final = sh_final in HP_TYPES
    dhw_hp_final = dhw_final in HP_TYPES

    sh_hp_initial = sh_initial in HP_TYPES
    dhw_hp_initial = dhw_initial in HP_TYPES

    if not sh_hp_final and not dhw_hp_final:
        return 0

    if sh_hp_initial and dhw_hp_initial and sh_initial == dhw_initial:
        old_cap = p95(sh_arr + dhw_arr)
        old_hp = pick_hp(old_cap, hp_catalog[sh_initial.split("_")[1]])
        old_cost = old_hp["cost_eur"]
    elif sh_hp_initial and dhw_hp_initial:
        old_cap_sh = p95(sh_arr)
        old_cap_dhw = p95(dhw_arr)
        old_hp_sh = pick_hp(old_cap_sh, hp_catalog[sh_initial.split("_")[1]])
        old_hp_dhw = pick_hp(old_cap_dhw, hp_catalog[dhw_initial.split("_")[1]])
        old_cost = old_hp_sh["cost_eur"] + old_hp_dhw["cost_eur"]
    elif sh_hp_initial:
        old_cap = p95(sh_arr)
        old_hp = pick_hp(old_cap, hp_catalog[sh_initial.split("_")[1]])
        old_cost = old_hp["cost_eur"]
    elif dhw_hp_initial:
        old_cap = p95(dhw_arr)
        old_hp = pick_hp(old_cap, hp_catalog[dhw_initial.split("_")[1]])
        old_cost = old_hp["cost_eur"]
    else:
        old_cost = 0

    if sh_hp_final and dhw_hp_final and sh_final == dhw_final:
        new_cap = p95(sh_arr + dhw_arr)
        new_hp = pick_hp(new_cap, hp_catalog[sh_final.split("_")[1]])
        new_cost = new_hp["cost_eur"]
    elif sh_hp_final and dhw_hp_final:
        new_cap_sh = p95(sh_arr)
        new_cap_dhw = p95(dhw_arr)
        new_hp_sh = pick_hp(new_cap_sh, hp_catalog[sh_final.split("_")[1]])
        new_hp_dhw = pick_hp(new_cap_dhw, hp_catalog[dhw_final.split("_")[1]])
        new_cost = new_hp_sh["cost_eur"] + new_hp_dhw["cost_eur"]
    elif sh_hp_final:
        new_cap = p95(sh_arr)
        new_hp = pick_hp(new_cap, hp_catalog[sh_final.split("_")[1]])
        new_cost = new_hp["cost_eur"]
    elif dhw_hp_final:
        new_cap = p95(dhw_arr)
        new_hp = pick_hp(new_cap, hp_catalog[dhw_final.split("_")[1]])
        new_cost = new_hp["cost_eur"]
    else:
        new_cost = 0

    return max(new_cost - old_cost, 0) * factor


def dhn_cost(building, interv, connection_cost):
    meta = building["meta"]

    sh_before = interv.get("sh_source", (meta["SHSource"], None))[0]
    dhw_before = interv.get("dhw_source", (meta["DHWsource"], None))[0]

    sh_after = interv.get("sh_source", (None, None))[1]
    dhw_after = interv.get("dhw_source", (None, None))[1]

    sh_final = sh_after if sh_after is not None else meta["SHSource"]
    dhw_final = dhw_after if dhw_after is not None else meta["DHWsource"]

    before_has_dhn = (sh_before == "dhn") or (dhw_before == "dhn")
    after_has_dhn = (sh_final == "dhn") or (dhw_final == "dhn")

    if after_has_dhn and not before_has_dhn:
        return connection_cost

    return 0


def compute_peak_cost(bought, time_index, price_per_kw):
    df = pd.DataFrame({
        "bought": bought,
        "time": time_index
    })

    df = df.set_index("time")

    monthly_peak = df["bought"].resample("M").max()
    monthly_cost = monthly_peak * price_per_kw

    hourly_cost = (
        monthly_cost.reindex(df.index, method="ffill")
        / df.resample("M").size().reindex(df.index, method="ffill")
    )

    return hourly_cost.values


def fuel_cost_array(fuel_name, demand_wh, fuels):
    f = fuels[fuel_name]

    var_price = f["price"]["variable"]["value"]
    fixed = f["price"]["fixed"]["value"]

    fixed_per_step = fixed / 8760

    return (demand_wh / 1000.0) * var_price + fixed_per_step


def hourly_dhn_cost(operation, p, area):
    demand = operation["dhn_bought"] / 1_000_000

    if demand.sum() == 0:
        return demand * 0

    area_fee = area * p["area fee per m2"]
    fixed_heat = p["fixed heat price per MWh"]
    var_heat = p["variable heat price per MWh"]
    admin = p["admin fee yearly"]
    sub_fixed = p["subscription fixed yearly per unit"]
    sub_var = p["subscription variable price per MWh"]
    vat = p["VAT"]

    hourly = demand * (var_heat + sub_var)
    hourly += fixed_heat / 8760

    if hasattr(area_fee, "values"):
        area_fee_value = area_fee.values[0]
    else:
        area_fee_value = area_fee

    total = hourly + (area_fee_value + admin + sub_fixed) / 8760
    total *= (1 + vat)

    return total


def hourly_grid_cost(operation, prices, pricing):
    df = prices.copy()

    df["time"] = pd.to_datetime(df["time"], utc=True)
    df = df.sort_values("time")

    mask = (df["time"] >= "2025-01-01") & (df["time"] < "2026-01-01")
    df = df.loc[mask]

    time = df["time"]
    price = df["price"].astype(float).values

    bought = np.array(operation["electricity_bought"], dtype=float) / 1000.0
    sold = np.array(operation["electricity_sold"], dtype=float) / 1000.0

    n = min(len(price), len(bought), len(sold))

    price = price[:n]
    time = time.iloc[:n]
    bought = bought[:n]
    sold = sold[:n]

    buy_cfg = pricing["buy"]
    sell_cfg = pricing["sell"]

    energy_tax = float(buy_cfg["energy tax per kWh"])
    cert = float(buy_cfg["electricity certificate cost per kWh"])
    grid_var = float(buy_cfg["grid local distribution cost monthly per kWh usage"])
    vat = float(buy_cfg["VAT"])

    grid_comp = float(sell_cfg["grid compensation"])
    tax_credit = float(sell_cfg["tax credit"])

    peak_price = float(buy_cfg["grid local distribution cost monthly per kW peak"])
    fixed_monthly = float(buy_cfg["grid local distribution cost monthly fix"])

    buy_price = (price + energy_tax + cert + grid_var) * (1 + vat)
    sell_price = price + grid_comp + tax_credit

    hourly_cost = bought * buy_price - sold * sell_price

    months = time.dt.to_period("M")
    monthly_add = np.zeros(len(hourly_cost))

    for m in months.unique():
        idx = (months == m).values
        peak_kw = np.max(bought[idx])
        hours = np.sum(idx)
        monthly_cost = peak_kw * peak_price + fixed_monthly
        monthly_add[idx] = monthly_cost / hours

    return hourly_cost + monthly_add


def _load_json_or_dict(obj, name="input"):
    if isinstance(obj, dict):
        return obj

    if isinstance(obj, (str, Path)):
        with open(obj, "r", encoding="utf-8") as f:
            return json.load(f)

    raise TypeError(
        f"{name} must be either a dict or a path-like object, got {type(obj).__name__}"
    )


def build_dict_gen(
    building_info,
    intervention_dict,
    baseline_gdf_path,
    configuration,
    ee_measure_path,
    pv_type_path,
    hp_catalog_path,
    grid_pricing_path,
    dhn_pricing_path,
    fuels_path,
    spot_price_path
):
    baseline_gdf = load_buildings(baseline_gdf_path)

    ee_measure = _load_json_or_dict(ee_measure_path, "ee_measure_path")
    pv_installs = _load_json_or_dict(pv_type_path, "pv_type_path")
    hp_catalog = _load_json_or_dict(hp_catalog_path, "hp_catalog_path")
    pricing = _load_json_or_dict(grid_pricing_path, "grid_pricing_path")
    dhn_pricing = _load_json_or_dict(dhn_pricing_path, "dhn_pricing_path")
    fuels = _load_json_or_dict(fuels_path, "fuels_path")

    prices = pd.read_csv(spot_price_path)
    prices["price"] = prices["price"].astype(float)

    buildings_dict = {}

    for idx, building in building_info.items():
        EUR_to_SEK = 10.86

        idx = int(idx)

        meta = {}
        interventions = intervention_dict.get(idx, {})
        operation = {}
        costs = {}

        row = baseline_gdf.loc[baseline_gdf["id"] == idx]

        if row.empty:
            raise ValueError(f"Building id {idx} not found in baseline_gdf")

        meta["PV_type"] = row["PVType"].values[0]

        pv_production = building["pv_production"]
        electricity_need = building["hp_electricity"] + building["base"]["appliance_electricity"]

        operation["electricity_bought"] = np.maximum(electricity_need - pv_production, 0)
        operation["electricity_sold"] = np.maximum(pv_production - electricity_need, 0)
        operation["dhn_bought"] = building["thermal"]["dhn_demand"]

        isfuel = (
            any(x in building["meta"]["DHWsource"] for x in ["gas", "oil", "bio"])
            or any(x in building["meta"]["SHSource"] for x in ["gas", "oil", "bio"])
        )

        if isfuel:
            fuel = next(
                (
                    f
                    for val in [building["meta"]["DHWsource"], building["meta"]["SHSource"]]
                    if isinstance(val, str) and "_" in val and "hp" not in val
                    for f in [val.split("_", 1)[1]]
                    if f in {"gas", "oil", "bio"}
                ),
                None
            )

            operation["fuel_bought"] = building["thermal"][f"{fuel}_demand"]
            meta["fuel"] = fuel

        area = building["base"]["pv_available_area"]
        floors = row["Floors"]
        used_area = area * floors

        efficiency_measure_cost = 0
        pv_install_cost = 0
        hp_install_cost = 0
        dhn_connection_cost = 0

        envelope_type = row["Envelope"].values[0]
        building_type = envelope_type[:3]

        costs["capital cost"] = {}
        costs["operational cost"] = {}

        if "envelope" in interventions:
            efficiency_measure = interventions["envelope"]

            before_measure = efficiency_measure[0]
            after_measure = efficiency_measure[1]

            if before_measure in ["medium", "deep"]:
                before_measure = before_measure + " " + building_type

            if after_measure in ["medium", "deep"]:
                after_measure = after_measure + " " + building_type

            wall_area = building["base"]["opaque_exposed_area"] - building["base"]["pv_available_area"]
            roof_area = building["base"]["pv_available_area"]
            window_area = building["base"]["glazing_area"]

            wall_cost_fix = (
                ee_measure[after_measure]["wall cost constant"]
                - ee_measure[before_measure]["wall cost constant"]
            )

            roof_cost_fix = (
                ee_measure[after_measure]["roof cost constant"]
                - ee_measure[before_measure]["roof cost constant"]
            )

            window_cost_fix = (
                ee_measure[after_measure]["window cost constant"]
                - ee_measure[before_measure]["window cost constant"]
            )

            wall_cost_persqm = (
                ee_measure[after_measure]["wall cost per square meter"]
                - ee_measure[before_measure]["wall cost per square meter"]
            )

            roof_cost_persqm = (
                ee_measure[after_measure]["roof cost per square meter"]
                - ee_measure[before_measure]["roof cost per square meter"]
            )

            window_cost_persqm = (
                ee_measure[after_measure]["window cost per square meter"]
                - ee_measure[before_measure]["window cost per square meter"]
            )

            efficiency_measure_cost = (
                wall_cost_fix
                + roof_cost_fix
                + window_cost_fix
                + wall_cost_persqm * wall_area
                + roof_cost_persqm * roof_area
                + window_cost_persqm * window_area
            )

        if "PVpercentage" in interventions:
            pv_install = interventions["PVpercentage"]

            before_measure = pv_install[0]
            after_measure = pv_install[1]

            pv_type = meta["PV_type"]

            pv_install_area = (
                (after_measure - before_measure)
                * building["base"]["pv_available_area"]
                / 100
            )

            pv_install_fix_cost = pv_installs[pv_type]["cost_fixed"]
            pv_install_persqm_cost = pv_installs[pv_type]["cost_per_m2"]

            if pv_install_area > 0:
                pv_install_cost = (
                    pv_install_persqm_cost * pv_install_area
                    + pv_install_fix_cost
                )

        if any(x in interventions for x in ["dhw_source", "sh_source"]):
            hp_install_cost = hp_cost(
                building,
                interventions,
                hp_catalog,
                factor=1.0
            ) * EUR_to_SEK

            dhn_connection_cost = dhn_cost(
                building,
                interventions,
                connection_cost=24000
            )

        costs["capital cost"]["efficiency_measure_cost"] = efficiency_measure_cost
        costs["capital cost"]["pv_install_cost"] = pv_install_cost
        costs["capital cost"]["heating_systems"] = hp_install_cost + dhn_connection_cost
        costs["capital cost"]["total"] = sum(costs["capital cost"].values())

        electricity_cost = hourly_grid_cost(operation, prices, pricing["1"])

        if len(electricity_cost) < 8760:
            electricity_cost = np.pad(
                electricity_cost,
                (8760 - len(electricity_cost), 0),
                mode="edge"
            )

        costs["operational cost"]["electricity"] = electricity_cost

        costs["operational cost"]["district_heating"] = hourly_dhn_cost(
            operation,
            dhn_pricing["1"]["buy"],
            used_area
        )

        if isfuel:
            costs["operational cost"]["fuel"] = (
                fuel_cost_array(meta["fuel"], operation["fuel_bought"], fuels)
                * EUR_to_SEK
            )

        costs["operational cost"]["hourly_total"] = np.sum(
            list(costs["operational cost"].values()),
            axis=0
        )

        costs["operational cost"]["yearly_total"] = np.sum(
            costs["operational cost"]["hourly_total"]
        )

        buildings_dict[idx] = {
            "meta": meta,
            "interventions": interventions,
            "operation": operation,
            "costs": costs
        }

    return buildings_dict