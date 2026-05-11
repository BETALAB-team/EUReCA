import os
import json
import numpy as np
from copy import deepcopy
from itertools import product


def copy_buildings(buildings):
    return deepcopy(buildings)


def copy_grid_levers(grid_levers):
    return deepcopy(grid_levers)


def copy_dhn_levers(dhn_levers):
    return deepcopy(dhn_levers)


def prepare_global_assumptions(assumptions):
    if assumptions is None:
        assumptions = {}

    tech_data = dict(assumptions)

    defaults = {
        "discount_rate": 0.04,
        "horizon_years": 20,
        "hp_scop": 3.0,
        "boiler_efficiency": 0.9,
        "pv_yield_kwh_per_m2_year": 160.0,
        "pv_self_consumption_ratio": 0.6,
        "pv_capex_per_m2": 2500.0,
        "hp_capex_sh": 90000.0,
        "hp_capex_dhw": 40000.0,
        "boiler_capex_sh": 35000.0,
        "boiler_capex_dhw": 20000.0,
        "generator_capex": 50000.0,
        "boiler_fuel_cost_per_kwh": 1.2,
        "generator_cost_per_kwh": 2.5,
        "generator_kwh_per_year": 0.0,
        "electricity_spot_price_per_kwh": 1.0,
        "dhn_connection_cost": 25000.0,
    }

    for k, v in defaults.items():
        if k not in tech_data:
            tech_data[k] = v

    return tech_data


def annualize_capex(capex, discount_rate, horizon_years):
    capex = float(capex)
    r = float(discount_rate)
    T = int(horizon_years)

    if capex == 0:
        return 0.0

    if r == 0:
        return capex / T

    crf = r * (1 + r) ** T / ((1 + r) ** T - 1)
    return capex * crf


def scale_grid_capex_to_regulated_area(
    study_area_grid_capex,
    study_area_customer_count,
    grid_regulation,
):
    study_area_grid_capex = float(study_area_grid_capex)
    study_area_customer_count = float(study_area_customer_count)

    if study_area_grid_capex <= 0:
        return 0.0

    if study_area_customer_count <= 0:
        return 0.0

    similarity = float(grid_regulation["background_similarity_factor"])
    background_customers = float(grid_regulation["background_baseline_customer_count"])

    represented_customers = study_area_customer_count + similarity * background_customers
    scale_factor = represented_customers / study_area_customer_count

    return study_area_grid_capex * scale_factor


def build_grid_map(grids):
    grid_map = {}

    for grid in grids:
        grid_id = grid.id

        if grid_id in grid_map:
            raise ValueError(f"Duplicate grid id found: {grid_id}")

        grid_map[grid_id] = grid

    return grid_map


def normalize_source(source):
    if source is None:
        return None

    s = str(source).lower()

    if "dhn" in s:
        return "dhn"
    if "hp" in s or "split" in s:
        return "hp"
    if "boiler" in s:
        return "boiler"

    return s


def summarize_building_base(building, tech_data):
    base = building["base"]

    sh = np.asarray(base["space_heating"], dtype=float)
    dhw = np.asarray(base["dhw_demand"], dtype=float)

    if "appliance_electricity" in base:
        elec = np.asarray(base["appliance_electricity"], dtype=float)
    elif "appliances" in base:
        elec = np.asarray(base["appliances"], dtype=float)
    elif "electricity" in base:
        elec = np.asarray(base["electricity"], dtype=float)
    else:
        raise KeyError("No electricity profile found in building['base']")

    pv_available_area = float(base.get("pv_available_area", 0.0))

    if "area" in base:
        area = float(base["area"])
    elif "net_floor_area" in base:
        area = float(base["net_floor_area"])
    elif "floor_area" in base:
        area = float(base["floor_area"])
    else:
        area = pv_available_area

    return {
        "space_heating_profile_wh": sh,
        "dhw_profile_wh": dhw,
        "electricity_profile_wh": elec,
        "space_heating_kwh_year": sh.sum() / 1000.0,
        "dhw_kwh_year": dhw.sum() / 1000.0,
        "electricity_kwh_year": elec.sum() / 1000.0,
        "space_heating_peak_kw": sh.max() / 1000.0,
        "dhw_peak_kw": dhw.max() / 1000.0,
        "electricity_peak_kw": elec.max() / 1000.0,
        "pv_available_area_m2": pv_available_area,
        "area_m2": area,
    }


def enumerate_building_choices(building, grid_map, dhns):
    reachable_grids = building.get("reachable_grids", [])
    reachable_dhns = building.get("reachable_dhns", [])

    reachable_grids = [g_id for g_id in reachable_grids if g_id in grid_map]
    reachable_dhns = [d_id for d_id in reachable_dhns if d_id in dhns]

    grid_choices = list(reachable_grids)
    dhn_choices = [None] + list(reachable_dhns)

    choice_space = []

    for sh_source, dhw_source, pv_installed, generator_installed, grid_id, dhn_id in product(
        ["dhn", "hp", "boiler"],
        ["dhn", "hp", "boiler"],
        [False, True],
        [False, True],
        grid_choices,
        dhn_choices,
    ):
        if sh_source == "dhn" or dhw_source == "dhn":
            if dhn_id is None:
                continue
        else:
            if dhn_id is not None:
                continue

        choice_space.append({
            "sh_source": sh_source,
            "dhw_source": dhw_source,
            "pv_installed": pv_installed,
            "generator_installed": generator_installed,
            "grid_id": grid_id,
            "dhn_id": dhn_id,
        })

    return choice_space


def calculate_choice_flows(base_summary, choice, tech_data):
    sh_kwh = float(base_summary["space_heating_kwh_year"])
    dhw_kwh = float(base_summary["dhw_kwh_year"])
    base_elec_kwh = float(base_summary["electricity_kwh_year"])
    pv_area_m2 = float(base_summary["pv_available_area_m2"])
    area_m2 = float(base_summary["area_m2"])
    electricity_peak_kw = float(base_summary["electricity_peak_kw"])

    hp_scop = float(tech_data["hp_scop"])
    pv_yield = float(tech_data["pv_yield_kwh_per_m2_year"])
    pv_self_consumption_ratio = float(tech_data["pv_self_consumption_ratio"])
    boiler_efficiency = float(tech_data["boiler_efficiency"])
    generator_kwh = float(tech_data["generator_kwh_per_year"])

    dhn_heat_kwh = 0.0
    hp_heat_kwh = 0.0
    boiler_heat_kwh = 0.0

    if choice["sh_source"] == "dhn":
        dhn_heat_kwh += sh_kwh
    elif choice["sh_source"] == "hp":
        hp_heat_kwh += sh_kwh
    elif choice["sh_source"] == "boiler":
        boiler_heat_kwh += sh_kwh
    else:
        raise ValueError(f"Unknown SH source: {choice['sh_source']}")

    if choice["dhw_source"] == "dhn":
        dhn_heat_kwh += dhw_kwh
    elif choice["dhw_source"] == "hp":
        hp_heat_kwh += dhw_kwh
    elif choice["dhw_source"] == "boiler":
        boiler_heat_kwh += dhw_kwh
    else:
        raise ValueError(f"Unknown DHW source: {choice['dhw_source']}")

    hp_electricity_kwh = hp_heat_kwh / hp_scop if hp_heat_kwh > 0 else 0.0
    boiler_fuel_kwh = boiler_heat_kwh / boiler_efficiency if boiler_heat_kwh > 0 else 0.0

    electricity_demand_kwh = base_elec_kwh + hp_electricity_kwh

    pv_generation_kwh = pv_area_m2 * pv_yield if choice["pv_installed"] else 0.0
    pv_used_on_site_kwh = min(electricity_demand_kwh, pv_generation_kwh * pv_self_consumption_ratio)
    pv_export_kwh = max(0.0, pv_generation_kwh - pv_used_on_site_kwh)

    generator_generation_kwh = generator_kwh if choice["generator_installed"] else 0.0

    electricity_bought_kwh = max(
        0.0,
        electricity_demand_kwh - pv_used_on_site_kwh - generator_generation_kwh
    )

    electricity_sold_kwh = max(
        0.0,
        pv_export_kwh + generator_generation_kwh + pv_used_on_site_kwh - electricity_demand_kwh
    )

    return {
        "area_m2": area_m2,
        "electricity_peak_kw_year": electricity_peak_kw,
        "space_heating_kwh_year": sh_kwh,
        "dhw_kwh_year": dhw_kwh,
        "base_electricity_kwh_year": base_elec_kwh,
        "dhn_heat_kwh_year": dhn_heat_kwh,
        "hp_heat_kwh_year": hp_heat_kwh,
        "boiler_heat_kwh_year": boiler_heat_kwh,
        "hp_electricity_kwh_year": hp_electricity_kwh,
        "boiler_fuel_kwh_year": boiler_fuel_kwh,
        "electricity_demand_kwh_year": electricity_demand_kwh,
        "pv_generation_kwh_year": pv_generation_kwh,
        "generator_generation_kwh_year": generator_generation_kwh,
        "electricity_bought_kwh_year": electricity_bought_kwh,
        "electricity_sold_kwh_year": electricity_sold_kwh,
    }


def extract_current_state(building):
    meta = building.get("meta", {})

    return {
        "current_grid_id": building.get("current_grid", building.get("current_grid_id")),
        "current_dhn_id": building.get("current_dhn", building.get("current_dhn_id")),
        "current_sh_source": normalize_source(meta.get("SHSource")),
        "current_dhw_source": normalize_source(meta.get("DHWsource")),
    }


def calculate_choice_capex_from_current_state(
    building,
    choice,
    current_state,
    tech_data,
):
    base = building["base"]
    pv_area_m2 = float(base.get("pv_available_area", 0.0))

    capex = {
        "hp_sh": 0.0,
        "hp_dhw": 0.0,
        "boiler_sh": 0.0,
        "boiler_dhw": 0.0,
        "pv": 0.0,
        "generator": 0.0,
        "dhn_connection": 0.0,
        "dhn_switch": 0.0,
        "grid_connection": 0.0,
        "grid_switch": 0.0,
    }

    current_sh_source = current_state["current_sh_source"]
    current_dhw_source = current_state["current_dhw_source"]
    current_dhn_id = current_state["current_dhn_id"]

    if choice["sh_source"] == "hp" and current_sh_source != "hp":
        capex["hp_sh"] = float(tech_data["hp_capex_sh"])
    elif choice["sh_source"] == "boiler" and current_sh_source != "boiler":
        capex["boiler_sh"] = float(tech_data["boiler_capex_sh"])

    if choice["dhw_source"] == "hp" and current_dhw_source != "hp":
        capex["hp_dhw"] = float(tech_data["hp_capex_dhw"])
    elif choice["dhw_source"] == "boiler" and current_dhw_source != "boiler":
        capex["boiler_dhw"] = float(tech_data["boiler_capex_dhw"])

    if choice["pv_installed"]:
        capex["pv"] = pv_area_m2 * float(tech_data["pv_capex_per_m2"])

    if choice["generator_installed"]:
        capex["generator"] = float(tech_data["generator_capex"])

    chosen_dhn_id = choice["dhn_id"]
    if chosen_dhn_id is not None and current_dhn_id is None:
        capex["dhn_connection"] = float(tech_data["dhn_connection_cost"])

    capex["total"] = sum(capex.values())

    return capex


def extract_customer_cost_factors(flows, choice):
    return {
        "dhn": {
            "active": choice["dhn_id"] is not None and float(flows["dhn_heat_kwh_year"]) > 0,
            "dhn_id": choice["dhn_id"],
            "area_m2": float(flows["area_m2"]),
            "heat_demand_mwh_year": float(flows["dhn_heat_kwh_year"]) / 1000.0,
        },
        "grid": {
            "active": choice["grid_id"] is not None,
            "grid_id": choice["grid_id"],
            "electricity_bought_kwh_year": float(flows["electricity_bought_kwh_year"]),
            "electricity_sold_kwh_year": float(flows["electricity_sold_kwh_year"]),
            "electricity_peak_kw_year": float(flows["electricity_peak_kw_year"]),
        },
        "boiler": {
            "active": float(flows["boiler_fuel_kwh_year"]) > 0,
            "boiler_fuel_kwh_year": float(flows["boiler_fuel_kwh_year"]),
        },
        "generator": {
            "active": choice["generator_installed"] and float(flows["generator_generation_kwh_year"]) > 0,
            "generator_generation_kwh_year": float(flows["generator_generation_kwh_year"]),
        },
    }


def extract_supplier_cost_factors(flows, choice, tech_data):
    return {
        "dhn": {
            "active": choice["dhn_id"] is not None and float(flows["dhn_heat_kwh_year"]) > 0,
            "dhn_id": choice["dhn_id"],
            "heat_sold_mwh_year": float(flows["dhn_heat_kwh_year"]) / 1000.0,
        },
        "grid": {
            "active": choice["grid_id"] is not None,
            "grid_id": choice["grid_id"],
            "electricity_bought_kwh_year": float(flows["electricity_bought_kwh_year"]),
            "electricity_sold_kwh_year": float(flows["electricity_sold_kwh_year"]),
        },
        "boiler": {
            "active": float(flows["boiler_fuel_kwh_year"]) > 0,
            "boiler_fuel_kwh_year": float(flows["boiler_fuel_kwh_year"]),
            "boiler_fuel_cost_per_kwh": float(tech_data["boiler_fuel_cost_per_kwh"]),
        },
        "generator": {
            "active": choice["generator_installed"] and float(flows["generator_generation_kwh_year"]) > 0,
            "generator_generation_kwh_year": float(flows["generator_generation_kwh_year"]),
            "generator_cost_per_kwh": float(tech_data["generator_cost_per_kwh"]),
        },
    }


def precompute_building_option(
    building,
    base_summary,
    choice,
    current_state,
    tech_data,
):
    flows = calculate_choice_flows(
        base_summary=base_summary,
        choice=choice,
        tech_data=tech_data,
    )

    capex = calculate_choice_capex_from_current_state(
        building=building,
        choice=choice,
        current_state=current_state,
        tech_data=tech_data,
    )

    customer_cost_factors = extract_customer_cost_factors(
        flows=flows,
        choice=choice,
    )

    supplier_cost_factors = extract_supplier_cost_factors(
        flows=flows,
        choice=choice,
        tech_data=tech_data,
    )

    return {
        "choice": dict(choice),
        "flows": flows,
        "capex": capex,
        "customer_cost_factors": customer_cost_factors,
        "supplier_cost_factors": supplier_cost_factors,
    }


def build_building_option_cache(
    buildings,
    grids,
    dhns,
    grid_levers,
    dhn_levers,
    assumptions,
):
    grid_map = build_grid_map(grids)
    tech_data = prepare_global_assumptions(assumptions)

    cache = {}

    for b_id, building in buildings.items():
        base_summary = summarize_building_base(building, tech_data)
        choice_space = enumerate_building_choices(building, grid_map, dhns)
        current_state = extract_current_state(building)

        options = []

        for choice in choice_space:
            option = precompute_building_option(
                building=building,
                base_summary=base_summary,
                choice=choice,
                current_state=current_state,
                tech_data=tech_data,
            )
            options.append(option)

        cache[b_id] = {
            "base_summary": base_summary,
            "options": options,
        }

    return cache


def calculate_choice_objective(capex, opex, tech_data):
    discount_rate = float(tech_data["discount_rate"])
    horizon_years = int(tech_data["horizon_years"])

    if discount_rate == 0:
        annuity_factor = horizon_years
    else:
        annuity_factor = (1 - (1 + discount_rate) ** (-horizon_years)) / discount_rate

    return {
        "capex_total": float(capex["total"]),
        "opex_yearly": float(opex["total"]),
        "npv_opex": float(opex["total"]) * annuity_factor,
        "npv_total": float(capex["total"]) + float(opex["total"]) * annuity_factor,
    }


def evaluate_cached_option(
    option,
    grid_levers,
    dhn_levers,
    tech_data,
    supplier_costs,
):
    choice = option["choice"]
    capex = option["capex"]
    customer_factors = option["customer_cost_factors"]

    customer_opex = {
        "dhn": 0.0,
        "grid": 0.0,
        "boiler": 0.0,
        "generator": 0.0,
        "total": 0.0,
    }

    supplier_contributions = {
        "grid": {
            "revenue": 0.0,
            "cost": 0.0,
            "profit": 0.0,
            "electricity_bought_kwh_year": 0.0,
            "electricity_sold_kwh_year": 0.0,
            "peak_kw": 0.0,
        },
        "dhn": {
            "revenue": 0.0,
            "cost": 0.0,
            "profit": 0.0,
            "heat_sold_mwh_year": 0.0,
        },
    }

    if customer_factors["dhn"]["active"]:
        dhn_id = customer_factors["dhn"]["dhn_id"]
        p = dhn_levers[dhn_id]["pricing"]

        area_m2 = float(customer_factors["dhn"]["area_m2"])
        demand_mwh = float(customer_factors["dhn"]["heat_demand_mwh_year"])

        area_fee = area_m2 * float(p["area fee per m2"])
        fixed_heat = float(p["fixed heat price per MWh"])
        var_heat = float(p["variable heat price per MWh"])
        admin = float(p["admin fee yearly"])
        sub_fixed = float(p["subscription fixed yearly per unit"])
        sub_var = float(p["subscription variable price per MWh"])
        vat = float(p["VAT"])

        dhn_customer_cost = demand_mwh * (var_heat + sub_var)
        dhn_customer_cost += fixed_heat
        dhn_customer_cost += area_fee + admin + sub_fixed
        dhn_customer_cost *= (1.0 + vat)

        customer_opex["dhn"] = dhn_customer_cost

        heat_supply_cost_per_mwh = float(
            supplier_costs["dhns"].get(dhn_id, {}).get("heat_supply_cost_per_mwh", 0.0)
        )

        supplier_contributions["dhn"]["revenue"] = dhn_customer_cost
        supplier_contributions["dhn"]["cost"] = demand_mwh * heat_supply_cost_per_mwh
        supplier_contributions["dhn"]["heat_sold_mwh_year"] = demand_mwh
        supplier_contributions["dhn"]["profit"] = (
            supplier_contributions["dhn"]["revenue"] -
            supplier_contributions["dhn"]["cost"]
        )

    if customer_factors["grid"]["active"]:
        grid_id = customer_factors["grid"]["grid_id"]
        p = grid_levers[grid_id]["pricing"]

        bought = float(customer_factors["grid"]["electricity_bought_kwh_year"])
        sold = float(customer_factors["grid"]["electricity_sold_kwh_year"])
        peak_kw = float(customer_factors["grid"]["electricity_peak_kw_year"])

        buy_cfg = p["buy"]
        sell_cfg = p["sell"]

        spot = float(tech_data["electricity_spot_price_per_kwh"])
        energy_tax = float(buy_cfg["energy tax per kWh"])
        cert = float(buy_cfg["electricity certificate cost per kWh"])
        grid_var = float(buy_cfg["grid local distribution cost monthly per kWh usage"])
        vat = float(buy_cfg["VAT"])

        grid_comp = float(sell_cfg["grid compensation"])
        tax_credit = float(sell_cfg["tax credit"])

        peak_price_monthly = float(buy_cfg["grid local distribution cost monthly per kW peak"])
        fixed_monthly = float(buy_cfg["grid local distribution cost monthly fix"])

        grid_customer_cost = bought * (spot + energy_tax + cert)
        grid_customer_cost += bought * grid_var
        grid_customer_cost += 12.0 * (peak_kw * peak_price_monthly + fixed_monthly)
        grid_customer_cost -= sold * (spot + grid_comp + tax_credit)
        grid_customer_cost *= (1.0 + vat)

        customer_opex["grid"] = grid_customer_cost

        electricity_purchase_cost_per_kwh = float(
            supplier_costs["grids"].get(grid_id, {}).get("electricity_purchase_cost_per_kwh", 0.0)
        )

        grid_fixed_cost_yearly = float(
            supplier_costs["grids"].get(grid_id, {}).get("grid_fixed_cost_yearly", 0.0)
        )

        supplier_contributions["grid"]["revenue"] = grid_customer_cost
        supplier_contributions["grid"]["cost"] = (
            bought * electricity_purchase_cost_per_kwh + grid_fixed_cost_yearly
        )
        supplier_contributions["grid"]["electricity_bought_kwh_year"] = bought
        supplier_contributions["grid"]["electricity_sold_kwh_year"] = sold
        supplier_contributions["grid"]["peak_kw"] = peak_kw
        supplier_contributions["grid"]["profit"] = (
            supplier_contributions["grid"]["revenue"] -
            supplier_contributions["grid"]["cost"]
        )

    if customer_factors["boiler"]["active"]:
        customer_opex["boiler"] = (
            float(customer_factors["boiler"]["boiler_fuel_kwh_year"]) *
            float(tech_data["boiler_fuel_cost_per_kwh"])
        )

    if customer_factors["generator"]["active"]:
        customer_opex["generator"] = (
            float(customer_factors["generator"]["generator_generation_kwh_year"]) *
            float(tech_data["generator_cost_per_kwh"])
        )

    customer_opex["total"] = (
        customer_opex["dhn"] +
        customer_opex["grid"] +
        customer_opex["boiler"] +
        customer_opex["generator"]
    )

    customer_objective = calculate_choice_objective(
        capex=capex,
        opex=customer_opex,
        tech_data=tech_data,
    )

    return {
        "choice": choice,
        "flows": option["flows"],
        "capex": capex,
        "customer_opex": customer_opex,
        "customer_objective": customer_objective,
        "supplier_contributions": supplier_contributions,
    }


def choose_best_cached_option_for_building(
    cached_building,
    grid_levers,
    dhn_levers,
    tech_data,
    supplier_costs,
):
    best_option_result = None
    best_objective = None

    for option in cached_building["options"]:
        option_result = evaluate_cached_option(
            option=option,
            grid_levers=grid_levers,
            dhn_levers=dhn_levers,
            tech_data=tech_data,
            supplier_costs=supplier_costs,
        )

        objective = float(option_result["customer_objective"]["npv_total"])

        if best_option_result is None or objective < best_objective:
            best_option_result = option_result
            best_objective = objective

    return best_option_result


def initialize_supplier_profit_state(grid_levers, dhn_levers):
    supplier_profits = {
        "grids": {},
        "dhns": {},
        "study_area_grid_totals": {
            "customer_count": 0.0,
            "electricity_bought_kwh_year": 0.0,
            "electricity_sold_kwh_year": 0.0,
            "peak_kw": 0.0,
        },
        "background_grid_totals": {
            "customer_count": 0.0,
            "electricity_bought_kwh_year": 0.0,
            "electricity_sold_kwh_year": 0.0,
            "peak_kw": 0.0,
        },
        "total_grid_totals": {
            "customer_count": 0.0,
            "electricity_bought_kwh_year": 0.0,
            "electricity_sold_kwh_year": 0.0,
            "peak_kw": 0.0,
            "revenue": 0.0,
            "cost": 0.0,
            "profit": 0.0,
        },
    }

    for grid_id in grid_levers:
        supplier_profits["grids"][grid_id] = {
            "customers": [],
            "revenue": 0.0,
            "cost": 0.0,
            "profit": 0.0,
            "electricity_bought_kwh_year": 0.0,
            "electricity_sold_kwh_year": 0.0,
            "peak_kw": 0.0,
            "customer_count": 0.0,
        }

    for dhn_id in dhn_levers:
        supplier_profits["dhns"][dhn_id] = {
            "customers": [],
            "revenue": 0.0,
            "cost": 0.0,
            "profit": 0.0,
            "heat_sold_mwh_year": 0.0,
        }

    return supplier_profits


def accumulate_option_into_supplier_profits(
    supplier_profits,
    building_id,
    option_result,
):
    choice = option_result["choice"]
    supplier_contrib = option_result["supplier_contributions"]

    grid_id = choice["grid_id"]

    if grid_id is not None:
        grid_entry = supplier_profits["grids"][grid_id]
        grid_entry["customers"].append(building_id)
        grid_entry["revenue"] += float(supplier_contrib["grid"]["revenue"])
        grid_entry["cost"] += float(supplier_contrib["grid"]["cost"])
        grid_entry["electricity_bought_kwh_year"] += float(
            supplier_contrib["grid"]["electricity_bought_kwh_year"]
        )
        grid_entry["electricity_sold_kwh_year"] += float(
            supplier_contrib["grid"]["electricity_sold_kwh_year"]
        )
        grid_entry["peak_kw"] += float(supplier_contrib["grid"]["peak_kw"])
        grid_entry["customer_count"] += 1.0

        supplier_profits["study_area_grid_totals"]["customer_count"] += 1.0
        supplier_profits["study_area_grid_totals"]["electricity_bought_kwh_year"] += float(
            supplier_contrib["grid"]["electricity_bought_kwh_year"]
        )
        supplier_profits["study_area_grid_totals"]["electricity_sold_kwh_year"] += float(
            supplier_contrib["grid"]["electricity_sold_kwh_year"]
        )
        supplier_profits["study_area_grid_totals"]["peak_kw"] += float(
            supplier_contrib["grid"]["peak_kw"]
        )

    dhn_id = choice["dhn_id"]

    if dhn_id is not None:
        dhn_entry = supplier_profits["dhns"][dhn_id]
        dhn_entry["customers"].append(building_id)
        dhn_entry["revenue"] += float(supplier_contrib["dhn"]["revenue"])
        dhn_entry["cost"] += float(supplier_contrib["dhn"]["cost"])
        dhn_entry["heat_sold_mwh_year"] += float(
            supplier_contrib["dhn"]["heat_sold_mwh_year"]
        )

    return supplier_profits


def calculate_background_grid_state_from_similarity(
    study_area_current_totals,
    study_area_baseline_grid_totals,
    grid_regulation,
):
    similarity = float(grid_regulation["background_similarity_factor"])

    background_baseline = {
        "customer_count": float(grid_regulation["background_baseline_customer_count"]),
        "electricity_bought_kwh_year": float(
            grid_regulation["background_baseline_electricity_bought_kwh_year"]
        ),
        "electricity_sold_kwh_year": float(
            grid_regulation["background_baseline_electricity_sold_kwh_year"]
        ),
        "peak_kw": float(grid_regulation["background_baseline_peak_kw"]),
    }

    background_current = {}

    for key in [
        "customer_count",
        "electricity_bought_kwh_year",
        "electricity_sold_kwh_year",
        "peak_kw",
    ]:
        baseline_study = float(study_area_baseline_grid_totals[key])
        current_study = float(study_area_current_totals[key])
        baseline_background = float(background_baseline[key])

        if baseline_study > 0:
            relative_change = (current_study - baseline_study) / baseline_study
        else:
            relative_change = 0.0

        background_current[key] = baseline_background * (1.0 + similarity * relative_change)

        if background_current[key] < 0.0:
            background_current[key] = 0.0

    return background_current


def finalize_supplier_profits(
    supplier_profits,
    supplier_costs,
    grid_regulation,
    study_area_baseline_grid_totals,
):
    study_area_current = supplier_profits["study_area_grid_totals"]

    background_current = calculate_background_grid_state_from_similarity(
        study_area_current_totals=study_area_current,
        study_area_baseline_grid_totals=study_area_baseline_grid_totals,
        grid_regulation=grid_regulation,
    )

    supplier_profits["background_grid_totals"] = background_current

    total_grid_totals = {
        "customer_count": (
            float(study_area_current["customer_count"]) +
            float(background_current["customer_count"])
        ),
        "electricity_bought_kwh_year": (
            float(study_area_current["electricity_bought_kwh_year"]) +
            float(background_current["electricity_bought_kwh_year"])
        ),
        "electricity_sold_kwh_year": (
            float(study_area_current["electricity_sold_kwh_year"]) +
            float(background_current["electricity_sold_kwh_year"])
        ),
        "peak_kw": (
            float(study_area_current["peak_kw"]) +
            float(background_current["peak_kw"])
        ),
        "revenue": 0.0,
        "cost": 0.0,
        "profit": 0.0,
    }

    for grid_id, entry in supplier_profits["grids"].items():
        entry["profit"] = float(entry["revenue"]) - float(entry["cost"])

    for dhn_id, entry in supplier_profits["dhns"].items():
        entry["profit"] = float(entry["revenue"]) - float(entry["cost"])

    if len(supplier_profits["grids"]) > 0:
        reference_grid_id = next(iter(supplier_profits["grids"]))
        reference_grid_entry = supplier_profits["grids"][reference_grid_id]

        study_area_customer_count = float(study_area_current["customer_count"])
        total_customer_count = float(total_grid_totals["customer_count"])

        if study_area_customer_count > 0:
            revenue_per_customer = float(reference_grid_entry["revenue"]) / study_area_customer_count
            cost_per_customer = float(reference_grid_entry["cost"]) / study_area_customer_count
        else:
            revenue_per_customer = 0.0
            cost_per_customer = 0.0

        annual_grid_capex_cost = float(grid_regulation.get("annual_grid_capex_cost", 0.0))

        total_grid_totals["revenue"] = revenue_per_customer * total_customer_count
        total_grid_totals["cost"] = (
            cost_per_customer * total_customer_count +
            annual_grid_capex_cost
        )
        total_grid_totals["profit"] = (
            total_grid_totals["revenue"] -
            total_grid_totals["cost"]
        )

    supplier_profits["total_grid_totals"] = total_grid_totals

    return supplier_profits


def evaluate_market_from_cache(
    building_option_cache,
    grid_levers,
    dhn_levers,
    tech_data,
    supplier_costs,
    grid_regulation,
    study_area_baseline_grid_totals,
):
    building_results = {}

    supplier_profits = initialize_supplier_profit_state(
        grid_levers=grid_levers,
        dhn_levers=dhn_levers,
    )

    for b_id, cached_building in building_option_cache.items():
        best_option_result = choose_best_cached_option_for_building(
            cached_building=cached_building,
            grid_levers=grid_levers,
            dhn_levers=dhn_levers,
            tech_data=tech_data,
            supplier_costs=supplier_costs,
        )

        building_results[b_id] = best_option_result

        supplier_profits = accumulate_option_into_supplier_profits(
            supplier_profits=supplier_profits,
            building_id=b_id,
            option_result=best_option_result,
        )

    supplier_profits = finalize_supplier_profits(
        supplier_profits=supplier_profits,
        supplier_costs=supplier_costs,
        grid_regulation=grid_regulation,
        study_area_baseline_grid_totals=study_area_baseline_grid_totals,
    )

    return {
        "building_results": building_results,
        "supplier_profits": supplier_profits,
    }


def calculate_study_area_grid_totals(building_results):
    totals = {
        "customer_count": 0.0,
        "electricity_bought_kwh_year": 0.0,
        "electricity_sold_kwh_year": 0.0,
        "peak_kw": 0.0,
    }

    for _, option_result in building_results.items():
        choice = option_result["choice"]
        supplier_contrib = option_result["supplier_contributions"]

        if choice["grid_id"] is None:
            continue

        totals["customer_count"] += 1.0
        totals["electricity_bought_kwh_year"] += float(
            supplier_contrib["grid"]["electricity_bought_kwh_year"]
        )
        totals["electricity_sold_kwh_year"] += float(
            supplier_contrib["grid"]["electricity_sold_kwh_year"]
        )
        totals["peak_kw"] += float(supplier_contrib["grid"]["peak_kw"])

    return totals


def list_market_actors(grid_levers, dhn_levers):
    actors = []

    if len(grid_levers) > 0:
        actors.append(("grid", "global"))

    if len(dhn_levers) > 0:
        actors.append(("dhn", "global"))

    return actors


def set_nested_value(d, dotted_key, value):
    keys = dotted_key.split(".")
    target = d

    for k in keys[:-1]:
        target = target[k]

    target[keys[-1]] = value


def get_nested_value(d, dotted_key):
    keys = dotted_key.split(".")
    value = d

    for k in keys:
        value = value[k]

    return value


def clip_to_bounds(value, bounds):
    lower, upper = bounds
    return max(lower, min(upper, value))


def generate_local_actor_candidates(
    actor_type,
    actor_id,
    grid_levers,
    dhn_levers,
    optimization_settings,
):
    local_search = optimization_settings["local_search"]

    if actor_type == "grid":
        if len(grid_levers) == 0:
            return []
        reference_id = next(iter(grid_levers))
        base_levers = deepcopy(grid_levers[reference_id])
    elif actor_type == "dhn":
        if len(dhn_levers) == 0:
            return []
        reference_id = next(iter(dhn_levers))
        base_levers = deepcopy(dhn_levers[reference_id])
    else:
        raise ValueError(f"Unknown actor_type: {actor_type}")

    actor_search = local_search[actor_type]
    active_levers = actor_search["active_levers"]
    step_sizes = actor_search["step_sizes"]
    bounds = actor_search["bounds"]

    candidates = [deepcopy(base_levers)]

    for lever_name in active_levers:
        current_value = float(get_nested_value(base_levers, lever_name))
        step = float(step_sizes[lever_name])
        lever_bounds = bounds[lever_name]

        for direction in (-1.0, 1.0):
            candidate = deepcopy(base_levers)
            new_value = clip_to_bounds(current_value + direction * step, lever_bounds)
            set_nested_value(candidate, lever_name, new_value)
            candidates.append(candidate)

    return candidates


def apply_best_response(
    grid_levers,
    dhn_levers,
    actor_type,
    actor_id,
    best_response,
):
    if isinstance(best_response, dict) and "levers" in best_response:
        levers = best_response["levers"]
    else:
        levers = best_response

    if actor_type == "grid":
        for grid_id in list(grid_levers.keys()):
            grid_levers[grid_id] = deepcopy(levers)

    elif actor_type == "dhn":
        for dhn_id in list(dhn_levers.keys()):
            dhn_levers[dhn_id] = deepcopy(levers)

    else:
        raise ValueError(f"Unknown actor_type: {actor_type}")

    return grid_levers, dhn_levers


def flatten_numeric(d, prefix=""):
    out = {}

    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else str(k)

        if isinstance(v, dict):
            out.update(flatten_numeric(v, key))
        elif isinstance(v, (int, float)):
            out[key] = float(v)

    return out


def compute_total_lever_change(
    old_grid_levers,
    old_dhn_levers,
    new_grid_levers,
    new_dhn_levers,
):
    total_change = 0.0

    all_grid_ids = set(old_grid_levers) | set(new_grid_levers)

    for grid_id in all_grid_ids:
        old_flat = flatten_numeric(old_grid_levers.get(grid_id, {}), prefix=f"grid.{grid_id}")
        new_flat = flatten_numeric(new_grid_levers.get(grid_id, {}), prefix=f"grid.{grid_id}")

        all_keys = set(old_flat) | set(new_flat)

        for key in all_keys:
            total_change += abs(new_flat.get(key, 0.0) - old_flat.get(key, 0.0))

    all_dhn_ids = set(old_dhn_levers) | set(new_dhn_levers)

    for dhn_id in all_dhn_ids:
        old_flat = flatten_numeric(old_dhn_levers.get(dhn_id, {}), prefix=f"dhn.{dhn_id}")
        new_flat = flatten_numeric(new_dhn_levers.get(dhn_id, {}), prefix=f"dhn.{dhn_id}")

        all_keys = set(old_flat) | set(new_flat)

        for key in all_keys:
            total_change += abs(new_flat.get(key, 0.0) - old_flat.get(key, 0.0))

    return total_change


def calculate_tariff_change_size(
    grid_levers,
    baseline_grid_levers,
    optimization_settings,
):
    if len(grid_levers) == 0 or len(baseline_grid_levers) == 0:
        return 0.0

    current_reference_id = next(iter(grid_levers))
    baseline_reference_id = next(iter(baseline_grid_levers))

    current_lever_set = grid_levers[current_reference_id]
    baseline_lever_set = baseline_grid_levers[baseline_reference_id]

    tariff_change_size = 0.0

    for lever_name in optimization_settings["local_search"]["grid"]["active_levers"]:
        current_value = float(get_nested_value(current_lever_set, lever_name))
        baseline_value = float(get_nested_value(baseline_lever_set, lever_name))
        tariff_change_size += (current_value - baseline_value) ** 2

    return tariff_change_size


def calculate_grid_regulation_score(
    market_state,
    grid_levers,
    baseline_grid_levers,
    grid_regulation,
    optimization_settings,
):
    total_grid_totals = market_state["supplier_profits"]["total_grid_totals"]

    base_allowed_revenue = float(grid_regulation["allowed_revenue_yearly"])
    annual_grid_capex_cost = float(grid_regulation.get("annual_grid_capex_cost", 0.0))

    allowed_revenue = base_allowed_revenue + annual_grid_capex_cost
    actual_revenue = float(total_grid_totals["revenue"])

    revenue_gap = actual_revenue - allowed_revenue

    tariff_change_size = calculate_tariff_change_size(
        grid_levers=grid_levers,
        baseline_grid_levers=baseline_grid_levers,
        optimization_settings=optimization_settings,
    )

    revenue_gap_weight = float(grid_regulation["revenue_gap_weight"])
    tariff_change_weight = float(grid_regulation["tariff_change_weight"])

    score = -(
        revenue_gap_weight * revenue_gap ** 2 +
        tariff_change_weight * tariff_change_size
    )

    return score


def optimize_single_actor_from_cache(
    actor_type,
    actor_id,
    building_option_cache,
    grid_levers,
    dhn_levers,
    tech_data,
    supplier_costs,
    optimization_settings,
    grid_regulation,
    baseline_grid_levers,
    study_area_baseline_grid_totals,
):
    candidates = generate_local_actor_candidates(
        actor_type=actor_type,
        actor_id=actor_id,
        grid_levers=grid_levers,
        dhn_levers=dhn_levers,
        optimization_settings=optimization_settings,
    )

    best_response = None
    best_score = None
    best_market_state = None

    for candidate in candidates:
        test_grid_levers = copy_grid_levers(grid_levers)
        test_dhn_levers = copy_dhn_levers(dhn_levers)

        test_grid_levers, test_dhn_levers = apply_best_response(
            grid_levers=test_grid_levers,
            dhn_levers=test_dhn_levers,
            actor_type=actor_type,
            actor_id=actor_id,
            best_response=candidate,
        )

        market_state = evaluate_market_from_cache(
            building_option_cache=building_option_cache,
            grid_levers=test_grid_levers,
            dhn_levers=test_dhn_levers,
            tech_data=tech_data,
            supplier_costs=supplier_costs,
            grid_regulation=grid_regulation,
            study_area_baseline_grid_totals=study_area_baseline_grid_totals,
        )

        if actor_type == "grid":
            score = calculate_grid_regulation_score(
                market_state=market_state,
                grid_levers=test_grid_levers,
                baseline_grid_levers=baseline_grid_levers,
                grid_regulation=grid_regulation,
                optimization_settings=optimization_settings,
            )

        elif actor_type == "dhn":
            score = sum(
                float(entry["profit"])
                for entry in market_state["supplier_profits"]["dhns"].values()
            )

        else:
            raise ValueError(f"Unknown actor_type: {actor_type}")

        if best_score is None or score > best_score:
            best_score = score
            best_response = candidate
            best_market_state = market_state

    return {
        "actor_type": actor_type,
        "actor_id": actor_id,
        "levers": best_response,
        "score": best_score,
        "market_state": best_market_state,
    }


def check_convergence(
    lever_change,
    history,
    optimization_settings,
):
    tolerance = float(optimization_settings.get("tolerance", 1e-6))
    min_rounds = int(optimization_settings.get("min_rounds", 1))

    if len(history) < min_rounds:
        return False

    return float(lever_change) <= tolerance


def optimize_market_levers(
    building_option_cache,
    initial_grid_levers,
    initial_dhn_levers,
    tech_data,
    supplier_costs,
    optimization_settings,
    grid_regulation,
    zero_study_area_baseline_grid_totals,
    avg_n_occ,
    study_area_grid_capex=0.0,
):
    grid_levers = copy_grid_levers(initial_grid_levers)
    dhn_levers = copy_dhn_levers(initial_dhn_levers)
    baseline_grid_levers = copy_grid_levers(initial_grid_levers)

    grid_regulation = deepcopy(grid_regulation)

    baseline_market_state = evaluate_market_from_cache(
        building_option_cache=building_option_cache,
        grid_levers=baseline_grid_levers,
        dhn_levers=dhn_levers,
        tech_data=tech_data,
        supplier_costs=supplier_costs,
        grid_regulation=grid_regulation,
        study_area_baseline_grid_totals=zero_study_area_baseline_grid_totals,
    )

    study_area_baseline_grid_totals = calculate_study_area_grid_totals(
        baseline_market_state["building_results"]
    )

    scaled_grid_capex = scale_grid_capex_to_regulated_area(
        study_area_grid_capex=study_area_grid_capex,
        study_area_customer_count=study_area_baseline_grid_totals["customer_count"],
        grid_regulation=grid_regulation,
    )

    annual_grid_capex_cost = annualize_capex(
        capex=scaled_grid_capex,
        discount_rate=tech_data["discount_rate"],
        horizon_years=tech_data["horizon_years"],
    )

    grid_regulation["study_area_grid_capex"] = float(study_area_grid_capex)
    grid_regulation["scaled_grid_capex"] = float(scaled_grid_capex)
    grid_regulation["annual_grid_capex_cost"] = float(annual_grid_capex_cost)

    history = []
    actors = list_market_actors(grid_levers, dhn_levers)

    for round_idx in range(optimization_settings["max_rounds"]):
        old_grid_levers = copy_grid_levers(grid_levers)
        old_dhn_levers = copy_dhn_levers(dhn_levers)

        round_history = {
            "round": round_idx,
            "actors": [],
        }

        for actor_type, actor_id in actors:
            best_response = optimize_single_actor_from_cache(
                actor_type=actor_type,
                actor_id=actor_id,
                building_option_cache=building_option_cache,
                grid_levers=grid_levers,
                dhn_levers=dhn_levers,
                tech_data=tech_data,
                supplier_costs=supplier_costs,
                optimization_settings=optimization_settings,
                grid_regulation=grid_regulation,
                baseline_grid_levers=baseline_grid_levers,
                study_area_baseline_grid_totals=study_area_baseline_grid_totals,
            )

            grid_levers, dhn_levers = apply_best_response(
                grid_levers=grid_levers,
                dhn_levers=dhn_levers,
                actor_type=actor_type,
                actor_id=actor_id,
                best_response=best_response,
            )

            round_history["actors"].append({
                "actor_type": actor_type,
                "actor_id": actor_id,
                "best_response": best_response,
            })

        market_state = evaluate_market_from_cache(
            building_option_cache=building_option_cache,
            grid_levers=grid_levers,
            dhn_levers=dhn_levers,
            tech_data=tech_data,
            supplier_costs=supplier_costs,
            grid_regulation=grid_regulation,
            study_area_baseline_grid_totals=study_area_baseline_grid_totals,
        )

        lever_change = compute_total_lever_change(
            old_grid_levers=old_grid_levers,
            old_dhn_levers=old_dhn_levers,
            new_grid_levers=grid_levers,
            new_dhn_levers=dhn_levers,
        )

        round_history["market_state"] = market_state
        round_history["lever_change"] = lever_change
        history.append(round_history)

        if check_convergence(
            lever_change=lever_change,
            history=history,
            optimization_settings=optimization_settings,
        ):
            break

    if avg_n_occ is not None and avg_n_occ != 0:
        for grid_id in grid_levers:
            grid_levers[grid_id]["pricing"]["buy"]["grid local distribution cost monthly fix"] /= avg_n_occ

    return {
        "grid_levers": grid_levers,
        "dhn_levers": dhn_levers,
        "history": history,
        "grid_regulation": grid_regulation,
        "study_area_baseline_grid_totals": study_area_baseline_grid_totals,
        "scaled_grid_capex": scaled_grid_capex,
        "annual_grid_capex_cost": annual_grid_capex_cost,
    }


def initialize_market_levers(
    dhns,
    grids,
    default_dhn_levers=None,
    default_grid_levers=None,
):
    if not isinstance(dhns, dict):
        raise TypeError("dhns must be a dictionary keyed by dhn_id")

    if not isinstance(grids, list):
        raise TypeError("grids must be a list of grid objects")

    if default_dhn_levers is None:
        default_dhn_levers = {
            "pricing": {
                "area fee per m2": 21.0,
                "fixed heat price per MWh": 140.0,
                "variable heat price per MWh": 384.0,
                "admin fee yearly": 500.0,
                "subscription fixed yearly per unit": 1415.0,
                "subscription variable price per MWh": 94.0,
                "VAT": 0.25,
            },
            "connection_cost": 25000.0,
            "financial": {
                "discount_rate": 0.04,
                "horizon_years": 20,
            },
        }

    if default_grid_levers is None:
        default_grid_levers = {
            "pricing": {
                "buy": {
                    "grid local distribution cost monthly fix": 100.0,
                    "grid local distribution cost monthly per kW peak": 78.0,
                    "grid local distribution cost monthly per kWh usage": 0.11,
                    "energy tax per kWh": 0.36,
                    "electricity certificate cost per kWh": 0.005,
                    "VAT": 0.25,
                },
                "sell": {
                    "grid compensation": 0.05,
                    "tax credit": 0.60,
                },
            },
            "financial": {
                "discount_rate": 0.04,
                "horizon_years": 20,
            },
        }

    dhn_levers = {}

    for dhn_id in dhns.keys():
        entry = deepcopy(default_dhn_levers)
        entry["dhn_id"] = dhn_id
        dhn_levers[dhn_id] = entry

    grid_levers = {}
    seen_grid_ids = set()

    for grid in grids:
        if not hasattr(grid, "id"):
            raise AttributeError("Each grid object must have an `id` attribute")

        grid_id = grid.id

        if grid_id in seen_grid_ids:
            raise ValueError(f"Duplicate grid id found: {grid_id}")

        seen_grid_ids.add(grid_id)

        entry = deepcopy(default_grid_levers)
        entry["grid_id"] = grid_id
        grid_levers[grid_id] = entry

    return dhn_levers, grid_levers


def build_network_maps(buildings, dhns, grids):
    for b_id in buildings:
        buildings[b_id]["current_grid"] = None
        buildings[b_id]["current_dhn"] = None

    for dhn_id, dhn in dhns.items():
        for node in dhn.nodes:
            if node.node_type == "consumer":
                b_id = str(node.node_id)

                if b_id in buildings:
                    buildings[b_id]["current_dhn"] = dhn_id

    for i, grid in enumerate(grids):
        grid.id = i
        grid_id = grid.id

        for node in grid.nodes:
            if node.node_type == "building":
                b_id = str(node.building_id)

                if b_id in buildings:
                    buildings[b_id]["current_grid"] = grid_id

    return buildings


def reachables(buildings, dhns, grids):
    for _, building in buildings.items():
        building["reachable_grids"] = [x for x, _ in enumerate(grids)]
        building["reachable_dhns"] = [x for x, _ in dhns.items()]

    return buildings