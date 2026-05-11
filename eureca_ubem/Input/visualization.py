import shutil
import webbrowser
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import folium

from shapely.geometry import LineString, Point
from folium.plugins import MeasureControl, Fullscreen


def array_total_kwh(x, unit="Wh"):
    if x is None:
        return 0.0

    arr = np.asarray(x, dtype=float)
    total = float(np.nansum(arr))

    if unit.lower() == "wh":
        return total / 1000.0

    if unit.lower() == "kwh":
        return total

    raise ValueError(f"unknown_operation_unit: {unit}")


def array_total_value(x):
    if x is None:
        return 0.0

    arr = np.asarray(x, dtype=float)
    return float(np.nansum(arr))


def safe_float(x, default=0.0):
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


def safe_divide(a, b, default=0.0):
    a = safe_float(a, default=0.0)
    b = safe_float(b, default=0.0)

    if b == 0:
        return default

    return a / b


def clean_fuel_name(fuel):
    if fuel is None:
        return ""

    fuel = str(fuel).strip().lower()

    replacements = {
        "biofuel": "bio",
        "biomass": "biomass",
        "natural gas": "natural_gas",
        "natgas": "natural_gas",
    }

    return replacements.get(fuel, fuel)


def get_nested(d, keys, default=None):
    current = d

    for key in keys:
        if not isinstance(current, dict):
            return default

        if key not in current:
            return default

        current = current[key]

    return current


def compute_total_dhn_generation(dhns, unit="kWh"):
    total_wh = 0.0

    if dhns is None:
        return 0.0

    iterator = dhns.values() if isinstance(dhns, dict) else dhns

    for dhn in iterator:
        hourly_heat_demand = getattr(dhn, "hourly_heat_demand", None)
        hourly_heat_loss = getattr(dhn, "hourly_heat_loss", None)

        if hourly_heat_demand is not None or hourly_heat_loss is not None:
            demand = 0.0 if hourly_heat_demand is None else np.asarray(hourly_heat_demand, dtype=float)
            loss = 0.0 if hourly_heat_loss is None else np.asarray(hourly_heat_loss, dtype=float)
            total_wh += float(np.nansum(demand + loss))
            continue

        nodes = getattr(dhn, "nodes", [])

        for node in nodes:
            if getattr(node, "node_type", None) != "supply":
                continue

            heat_injected = getattr(node, "heat_injected", None)

            if heat_injected is None:
                continue

            total_wh += float(np.nansum(np.asarray(heat_injected, dtype=float)))

    unit = unit.lower()

    if unit == "wh":
        return total_wh

    if unit == "kwh":
        return total_wh / 1_000.0

    if unit == "mwh":
        return total_wh / 1_000_000.0

    if unit == "gwh":
        return total_wh / 1_000_000_000.0

    raise ValueError(f"unknown_dhn_generation_unit: {unit}")

def is_active_line(line):
    active = getattr(line, "active", True)

    if isinstance(active, str):
        return active.strip().lower() in {"true", "1", "yes", "y"}

    return bool(active)

def compute_dhn_operation_from_objects(dhns, unit="Wh"):
    total_demand_kwh = 0.0
    total_loss_kwh = 0.0
    total_generation_kwh = 0.0

    if dhns is None:
        return {
            "dhn_heat_demand_kwh": 0.0,
            "dhn_heat_loss_kwh": 0.0,
            "dhn_heat_generation_kwh": 0.0,
        }

    iterator = dhns.values() if isinstance(dhns, dict) else dhns

    for dhn in iterator:
        demand = getattr(dhn, "hourly_heat_demand", None)
        loss = getattr(dhn, "hourly_heat_loss", None)

        demand_kwh = array_total_kwh(demand, unit=unit) if demand is not None else 0.0
        loss_kwh = array_total_kwh(loss, unit=unit) if loss is not None else 0.0

        total_demand_kwh += demand_kwh
        total_loss_kwh += loss_kwh
        total_generation_kwh += demand_kwh + loss_kwh

    return {
        "dhn_heat_demand_kwh": total_demand_kwh,
        "dhn_heat_loss_kwh": total_loss_kwh,
        "dhn_heat_generation_kwh": total_generation_kwh,
    }


def source_is_dhn(x):
    if x is None:
        return False

    x = str(x).strip().lower()

    return x in {
        "dhn",
        "district_heating",
        "district heating",
        "district-heating",
    }


def source_is_hp(x):
    if x is None:
        return False

    x = str(x).strip().lower()

    return (
        x in {
            "hp",
            "heatpump",
            "heat_pump",
            "heat pump",
            "hp_le",
            "hp_me",
            "hp_he",
            "le_heatpump",
            "me_heatpump",
            "he_heatpump",
            "le heatpump",
            "me heatpump",
            "he heatpump",
        }
        or x.startswith("hp_")
        or x.endswith("_heatpump")
    )


def count_building_heat_sources(config, buildings_gdf=None):
    dhn_connected = 0
    hp_connected = 0

    if isinstance(config, dict):
        iterator = config.values()

        for c in iterator:
            if not isinstance(c, dict):
                continue

            sh = (
                c.get("SHSource")
                or c.get("sh_source")
                or c.get("SHsource")
                or c.get("space_heating_source")
            )
            dhw = (
                c.get("DHWsource")
                or c.get("DHWSource")
                or c.get("dhw_source")
                or c.get("domestic_hot_water_source")
            )

            if source_is_dhn(sh) or source_is_dhn(dhw):
                dhn_connected += 1

            if source_is_hp(sh) or source_is_hp(dhw):
                hp_connected += 1

        return {
            "dhn_connected_buildings": dhn_connected,
            "hp_buildings": hp_connected,
        }

    if buildings_gdf is not None:
        for _, row in buildings_gdf.iterrows():
            sh = row.get("SHSource", None)
            dhw = row.get("DHWsource", None)

            if source_is_dhn(sh) or source_is_dhn(dhw):
                dhn_connected += 1

            if source_is_hp(sh) or source_is_hp(dhw):
                hp_connected += 1

    return {
        "dhn_connected_buildings": dhn_connected,
        "hp_buildings": hp_connected,
    }


def summarize_one_optimal_building(opt_data, assumptions):
    operation = opt_data.get("operation", {})
    meta = opt_data.get("meta", {})
    costs = opt_data.get("costs", {})

    unit = assumptions.get("operation_unit", "Wh")
    pef = assumptions.get("primary_energy_factors", {})

    electricity_bought_kwh = array_total_kwh(operation.get("electricity_bought"), unit)
    electricity_sold_kwh = array_total_kwh(operation.get("electricity_sold"), unit)
    dhn_bought_kwh = array_total_kwh(operation.get("dhn_bought"), unit)
    fuel_bought_kwh = array_total_kwh(operation.get("fuel_bought"), unit)

    fuel = clean_fuel_name(meta.get("fuel"))
    fuel_factor = pef.get(fuel, pef.get("fuel_default", 1.0))

    primary_energy_kwh = (
        electricity_bought_kwh * pef.get("electricity_bought", 1.0)
        + dhn_bought_kwh * pef.get("dhn_bought", 1.0)
        + fuel_bought_kwh * fuel_factor
    )

    capital_cost = 0.0
    npv = 0.0
    electricity_operational_cost = 0.0
    dhn_operational_cost = 0.0
    fuel_operational_cost = 0.0
    total_operational_cost = 0.0

    if isinstance(costs, dict):
        capital_cost_dict = costs.get("capital cost", {})
        financial_dict = costs.get("financial", {})
        operational_cost_dict = costs.get("operational cost", {})

        if isinstance(capital_cost_dict, dict):
            capital_cost = safe_float(capital_cost_dict.get("total", 0.0))

        if isinstance(financial_dict, dict):
            npv = safe_float(financial_dict.get("NPV", 0.0))

        if isinstance(operational_cost_dict, dict):
            electricity_operational_cost = array_total_value(
                operational_cost_dict.get("electricity", 0.0)
            )
            dhn_operational_cost = array_total_value(
                operational_cost_dict.get("district_heating", 0.0)
            )
            fuel_operational_cost = array_total_value(
                operational_cost_dict.get("fuel", 0.0)
            )

            if "yearly_total" in operational_cost_dict:
                total_operational_cost = safe_float(
                    operational_cost_dict.get("yearly_total", 0.0)
                )
            else:
                total_operational_cost = (
                    electricity_operational_cost
                    + dhn_operational_cost
                    + fuel_operational_cost
                )

    return {
        "electricity_bought_kwh": electricity_bought_kwh,
        "electricity_sold_kwh": electricity_sold_kwh,
        "dhn_bought_kwh": dhn_bought_kwh,
        "fuel_bought_kwh": fuel_bought_kwh,
        "fuel_type": fuel,
        "primary_energy_kwh": primary_energy_kwh,
        "building_capital_cost": capital_cost,
        "building_npv": npv,
        "electricity_operational_cost": electricity_operational_cost,
        "dhn_operational_cost": dhn_operational_cost,
        "fuel_operational_cost": fuel_operational_cost,
        "total_operational_cost": total_operational_cost,
    }


def summarize_optimal_set(optimal_set, assumptions):
    summaries = {}

    if optimal_set is None:
        return summaries

    for building_id, opt_data in optimal_set.items():
        if not isinstance(opt_data, dict):
            continue

        summary = summarize_one_optimal_building(opt_data, assumptions)

        summaries[building_id] = summary
        summaries[str(building_id)] = summary

        try:
            summaries[int(building_id)] = summary
        except Exception:
            pass

    return summaries


def compute_building_area_columns(gdf):
    out = gdf.copy()

    if out.crs is None:
        out = out.set_crs("EPSG:3006")

    if out.crs.is_geographic:
        area_gdf = out.to_crs("EPSG:3006")
    else:
        area_gdf = out

    out["footprint_area_m2"] = area_gdf.geometry.area.astype(float)

    if "Floors" in out.columns:
        floors = pd.to_numeric(out["Floors"], errors="coerce").fillna(1.0)
    else:
        floors = pd.Series(1.0, index=out.index)

    floors = floors.clip(lower=1.0)

    out["total_floor_area_m2"] = out["footprint_area_m2"] * floors

    return out


def attach_energy_attributes_to_buildings(buildings_gdf, optimal_set, assumptions):
    gdf = buildings_gdf.copy()
    gdf = compute_building_area_columns(gdf)

    summaries = summarize_optimal_set(optimal_set, assumptions)

    numeric_cols = [
        "electricity_bought_kwh",
        "electricity_sold_kwh",
        "dhn_bought_kwh",
        "fuel_bought_kwh",
        "primary_energy_kwh",
        "building_capital_cost",
        "building_npv",
        "electricity_operational_cost",
        "dhn_operational_cost",
        "fuel_operational_cost",
        "total_operational_cost",
        "npv_per_m2",
    ]

    for col in numeric_cols:
        if col not in gdf.columns:
            gdf[col] = 0.0

    if "fuel_type" not in gdf.columns:
        gdf["fuel_type"] = ""

    for idx, row in gdf.iterrows():
        candidates = []

        if "id" in gdf.columns:
            candidates.append(row["id"])
            candidates.append(str(row["id"]))

            try:
                candidates.append(int(row["id"]))
            except Exception:
                pass

        if "Name" in gdf.columns:
            candidates.append(row["Name"])
            candidates.append(str(row["Name"]))

        summary = None

        for c in candidates:
            if c in summaries:
                summary = summaries[c]
                break

        if summary is None:
            continue

        for k, v in summary.items():
            gdf.at[idx, k] = v

    area = pd.to_numeric(gdf["total_floor_area_m2"], errors="coerce").replace(0, np.nan)
    npv = pd.to_numeric(gdf["building_npv"], errors="coerce").fillna(0.0)

    gdf["npv_per_m2"] = (
        npv / area
    ).replace([np.inf, -np.inf], np.nan).fillna(0.0)

    return gdf


def summarize_map_totals(buildings_gdf):
    totals = {}

    for col in [
        "electricity_bought_kwh",
        "electricity_sold_kwh",
        "dhn_bought_kwh",
        "fuel_bought_kwh",
        "primary_energy_kwh",
        "building_capital_cost",
        "building_npv",
        "electricity_operational_cost",
        "dhn_operational_cost",
        "fuel_operational_cost",
        "total_operational_cost",
        "footprint_area_m2",
        "total_floor_area_m2",
    ]:
        if col in buildings_gdf.columns:
            totals[col] = float(
                pd.to_numeric(buildings_gdf[col], errors="coerce")
                .fillna(0.0)
                .sum()
            )
        else:
            totals[col] = 0.0

    if totals["total_floor_area_m2"] > 0:
        totals["npv_per_m2_average"] = (
            totals["building_npv"] / totals["total_floor_area_m2"]
        )
    else:
        totals["npv_per_m2_average"] = 0.0

    totals["average_grid_price_sek_per_kwh"] = safe_divide(
        totals["electricity_operational_cost"],
        totals["electricity_bought_kwh"],
    )

    totals["average_dhn_price_from_building_operation_sek_per_kwh"] = safe_divide(
        totals["dhn_operational_cost"],
        totals["dhn_bought_kwh"],
    )

    return totals


def safe_total_cost(x):
    if x is None:
        return 0.0

    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)

    if isinstance(x, dict):
        keys = [
            "total_cost",
            "total_cost_sek",
            "capex",
            "capex_sek",
            "cost",
            "cost_sek",
            "total",
        ]

        for k in keys:
            if k in x and isinstance(x[k], (int, float, np.integer, np.floating)):
                return float(x[k])

        total = 0.0

        for v in x.values():
            total += safe_total_cost(v)

        return total

    if isinstance(x, pd.DataFrame):
        numeric = x.select_dtypes(include=[np.number])

        if numeric.empty:
            return 0.0

        return float(numeric.sum().sum())

    return 0.0


def object_key(obj, names):
    if isinstance(obj, (str, int, float, np.integer, np.floating)):
        return obj

    for name in names:
        if hasattr(obj, name):
            return getattr(obj, name)

    return obj


def coords_look_like_lonlat(xs, ys):
    if len(xs) == 0 or len(ys) == 0:
        return False

    xs = np.array(xs, dtype=float)
    ys = np.array(ys, dtype=float)

    return (
        np.nanmin(xs) >= -180
        and np.nanmax(xs) <= 180
        and np.nanmin(ys) >= -90
        and np.nanmax(ys) <= 90
    )


def infer_xy_crs_from_buildings(buildings_gdf):
    if buildings_gdf is not None and buildings_gdf.crs is not None:
        return buildings_gdf.crs

    return "EPSG:3006"


def count_container(x):
    if x is None:
        return 0

    if isinstance(x, dict):
        return len(x)

    try:
        return len(x)
    except Exception:
        return 1


def dhn_to_gdf(retrofit, buildings_gdf):
    records = []
    xs = []
    ys = []

    dhn_systems = getattr(retrofit, "District_Heating_Systems", None)

    if dhn_systems is None:
        return None

    iterator = dhn_systems.items() if isinstance(dhn_systems, dict) else enumerate(dhn_systems)

    for dhn_id, dhn_obj in iterator:
        nodes = getattr(dhn_obj, "nodes", [])
        lines = getattr(dhn_obj, "lines", [])

        node_map = {}

        for node in nodes:
            node_id = getattr(node, "node_id", None)
            x = getattr(node, "x", None)
            y = getattr(node, "y", None)

            if node_id is None or x is None or y is None:
                continue

            x = float(x)
            y = float(y)

            node_map[node_id] = (x, y)
            node_map[str(node_id)] = (x, y)

            try:
                node_map[int(node_id)] = (x, y)
            except Exception:
                pass

            xs.append(x)
            ys.append(y)

        for line in lines:
            if not is_active_line(line):
                continue
            line_id = getattr(line, "line_id", None)

            start_node = object_key(getattr(line, "start_node", None), ["node_id"])
            end_node = object_key(getattr(line, "end_node", None), ["node_id"])

            if start_node not in node_map and str(start_node) in node_map:
                start_node = str(start_node)

            if end_node not in node_map and str(end_node) in node_map:
                end_node = str(end_node)

            if start_node not in node_map:
                try:
                    if int(start_node) in node_map:
                        start_node = int(start_node)
                except Exception:
                    pass

            if end_node not in node_map:
                try:
                    if int(end_node) in node_map:
                        end_node = int(end_node)
                except Exception:
                    pass

            if start_node not in node_map or end_node not in node_map:
                continue

            p1 = node_map[start_node]
            p2 = node_map[end_node]

            if p1 == p2:
                continue

            records.append(
                {
                    "network_type": "DHN",
                    "dhn_id": str(dhn_id),
                    "line_id": str(line_id),
                    "start_node": str(start_node),
                    "end_node": str(end_node),
                    "geometry": LineString([p1, p2]),
                }
            )

    if not records:
        return None

    crs = "EPSG:4326" if coords_look_like_lonlat(xs, ys) else infer_xy_crs_from_buildings(buildings_gdf)

    return gpd.GeoDataFrame(records, geometry="geometry", crs=crs)


def grid_to_gdf(retrofit, buildings_gdf):
    records = []
    xs = []
    ys = []

    grids = getattr(retrofit, "Electrical_Network", None)

    if grids is None:
        return None

    iterator = grids.items() if isinstance(grids, dict) else enumerate(grids)

    for grid_id, grid_obj in iterator:
        nodes = getattr(grid_obj, "nodes", [])
        lines = getattr(grid_obj, "lines", [])

        node_map = {}

        for node in nodes:
            bus_id = getattr(node, "bus_id", None)
            x = getattr(node, "x", None)
            y = getattr(node, "y", None)

            if bus_id is None or x is None or y is None:
                continue

            x = float(x)
            y = float(y)

            node_map[bus_id] = (x, y)
            node_map[str(bus_id)] = (x, y)

            try:
                node_map[int(bus_id)] = (x, y)
            except Exception:
                pass

            xs.append(x)
            ys.append(y)

        for line in lines:
            line_id = getattr(line, "line_id", None)

            from_bus = object_key(getattr(line, "from_bus", None), ["bus_id"])
            to_bus = object_key(getattr(line, "to_bus", None), ["bus_id"])

            if from_bus not in node_map and str(from_bus) in node_map:
                from_bus = str(from_bus)

            if to_bus not in node_map and str(to_bus) in node_map:
                to_bus = str(to_bus)

            if from_bus not in node_map:
                try:
                    if int(from_bus) in node_map:
                        from_bus = int(from_bus)
                except Exception:
                    pass

            if to_bus not in node_map:
                try:
                    if int(to_bus) in node_map:
                        to_bus = int(to_bus)
                except Exception:
                    pass

            if from_bus not in node_map or to_bus not in node_map:
                continue

            p1 = node_map[from_bus]
            p2 = node_map[to_bus]

            if p1 == p2:
                continue

            records.append(
                {
                    "network_type": "GRID",
                    "grid_id": str(grid_id),
                    "line_id": str(line_id),
                    "from_bus": str(from_bus),
                    "to_bus": str(to_bus),
                    "geometry": LineString([p1, p2]),
                }
            )

    if not records:
        return None

    crs = "EPSG:4326" if coords_look_like_lonlat(xs, ys) else infer_xy_crs_from_buildings(buildings_gdf)

    return gpd.GeoDataFrame(records, geometry="geometry", crs=crs)


def nodes_to_gdf(retrofit, buildings_gdf):
    records = []
    xs = []
    ys = []

    dhn_systems = getattr(retrofit, "District_Heating_Systems", None)

    if dhn_systems is not None:
        iterator = dhn_systems.items() if isinstance(dhn_systems, dict) else enumerate(dhn_systems)

        for dhn_id, dhn_obj in iterator:
            for node in getattr(dhn_obj, "nodes", []):
                node_id = getattr(node, "node_id", None)
                x = getattr(node, "x", None)
                y = getattr(node, "y", None)

                if node_id is None or x is None or y is None:
                    continue

                x = float(x)
                y = float(y)

                records.append(
                    {
                        "network_type": "DHN node",
                        "network_id": str(dhn_id),
                        "node_id": str(node_id),
                        "node_type": str(getattr(node, "node_type", "")),
                        "x": x,
                        "y": y,
                        "geometry": Point(x, y),
                    }
                )

                xs.append(x)
                ys.append(y)

    grids = getattr(retrofit, "Electrical_Network", None)

    if grids is not None:
        iterator = grids.items() if isinstance(grids, dict) else enumerate(grids)

        for grid_id, grid_obj in iterator:
            for node in getattr(grid_obj, "nodes", []):
                bus_id = getattr(node, "bus_id", None)
                x = getattr(node, "x", None)
                y = getattr(node, "y", None)

                if bus_id is None or x is None or y is None:
                    continue

                x = float(x)
                y = float(y)

                records.append(
                    {
                        "network_type": "Grid node",
                        "network_id": str(grid_id),
                        "node_id": str(bus_id),
                        "node_type": "",
                        "x": x,
                        "y": y,
                        "geometry": Point(x, y),
                    }
                )

                xs.append(x)
                ys.append(y)

    if not records:
        return None

    crs = "EPSG:4326" if coords_look_like_lonlat(xs, ys) else infer_xy_crs_from_buildings(buildings_gdf)

    return gpd.GeoDataFrame(records, geometry="geometry", crs=crs)


def as_wgs84(gdf):
    if gdf is None:
        return None

    if len(gdf) == 0:
        return None

    gdf = gdf.copy()

    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:3006")

    return gdf.to_crs("EPSG:4326")


def map_center_from_layers(*layers):
    valid = []

    for layer in layers:
        wgs = as_wgs84(layer)

        if wgs is not None and len(wgs) > 0:
            valid.append(wgs)

    if not valid:
        return [59.8586, 17.6389]

    merged = pd.concat(valid, ignore_index=True)
    merged = gpd.GeoDataFrame(merged, geometry="geometry", crs="EPSG:4326")

    centroid = merged.geometry.unary_union.centroid

    return [centroid.y, centroid.x]


def format_number(x, digits=2):
    try:
        return f"{float(x):,.{digits}f}"
    except Exception:
        return str(x)


def make_side_panel_html(
    game_step,
    capital_costs_grid,
    capital_costs_dhn,
    changes,
    dhn_count,
    dhn_line_count,
    grid_count,
    grid_line_count,
    building_count,
    energy_totals,
    dhn_operation_totals,
    source_counts,
):
    grid_capex_total = safe_total_cost(capital_costs_grid)
    dhn_capex_total = safe_total_cost(capital_costs_dhn)

    dhn_heat_demand_kwh = dhn_operation_totals.get("dhn_heat_demand_kwh", 0.0)

    avg_grid_price = safe_divide(
        energy_totals.get("electricity_operational_cost", 0.0),
        energy_totals.get("electricity_bought_kwh", 0.0),
    )

    avg_dhn_price = safe_divide(
        energy_totals.get("dhn_operational_cost", 0.0),
        dhn_heat_demand_kwh,
    )

    if changes:
        change_html = "<h3>Building configuration changes</h3><table>"

        for key_path, old_value, new_value in changes[:40]:
            change_html += f"""
            <tr>
                <td>{key_path}</td>
                <td>{old_value}</td>
                <td>→</td>
                <td>{new_value}</td>
            </tr>
            """

        change_html += "</table>"

        if len(changes) > 40:
            change_html += f"<p>Showing first 40 of {len(changes)} changes.</p>"
    else:
        change_html = "<h3>Building configuration changes</h3><p>No building changed configuration.</p>"

    html = f"""
    <div id="market-panel">
        <h2>Game step {game_step}</h2>

        <h3>Map objects</h3>
        <table>
            <tr><td>Buildings</td><td>{building_count}</td></tr>
            <tr><td>DHN systems</td><td>{dhn_count}</td></tr>
            <tr><td>DHN pipes</td><td>{dhn_line_count}</td></tr>
            <tr><td>Grid systems</td><td>{grid_count}</td></tr>
            <tr><td>Grid cable lines</td><td>{grid_line_count}</td></tr>
        </table>

        <h3>Average operational prices</h3>
        <table>
            <tr><td>Average grid price</td><td>{format_number(avg_grid_price, 4)} SEK/kWh</td></tr>
            <tr><td>Average DHN price</td><td>{format_number(avg_dhn_price, 4)} SEK/kWh</td></tr>
        </table>

        <h3>Total yearly operation</h3>
        <table>
            <tr><td>Grid electricity bought</td><td>{format_number(energy_totals.get("electricity_bought_kwh", 0.0), 2)} kWh</td></tr>
            <tr><td>Grid electricity sold</td><td>{format_number(energy_totals.get("electricity_sold_kwh", 0.0), 2)} kWh</td></tr>
            <tr><td>DHN heat demand</td><td>{format_number(dhn_operation_totals.get("dhn_heat_demand_kwh", 0.0), 2)} kWh</td></tr>
            <tr><td>DHN heat loss</td><td>{format_number(dhn_operation_totals.get("dhn_heat_loss_kwh", 0.0), 2)} kWh</td></tr>
            <tr><td>DHN heat generation</td><td>{format_number(dhn_operation_totals.get("dhn_heat_generation_kwh", 0.0), 2)} kWh</td></tr>
            <tr><td>Fuel bought</td><td>{format_number(energy_totals.get("fuel_bought_kwh", 0.0), 2)} kWh</td></tr>
            <tr><td>Primary energy</td><td>{format_number(energy_totals.get("primary_energy_kwh", 0.0), 2)} kWh PE</td></tr>
        </table>

        <h3>Total yearly operational costs</h3>
        <table>
            <tr><td>Electricity operational cost</td><td>{format_number(energy_totals.get("electricity_operational_cost", 0.0), 2)} SEK</td></tr>
            <tr><td>DHN operational cost</td><td>{format_number(energy_totals.get("dhn_operational_cost", 0.0), 2)} SEK</td></tr>
            <tr><td>Fuel operational cost</td><td>{format_number(energy_totals.get("fuel_operational_cost", 0.0), 2)} SEK</td></tr>
            <tr><td>Total operational cost</td><td>{format_number(energy_totals.get("total_operational_cost", 0.0), 2)} SEK</td></tr>
        </table>

        <h3>Building economics</h3>
        <table>
            <tr><td>Total building CAPEX</td><td>{format_number(energy_totals.get("building_capital_cost", 0.0), 2)} SEK</td></tr>
            <tr><td>Total building NPV</td><td>{format_number(energy_totals.get("building_npv", 0.0), 2)} SEK</td></tr>
            <tr><td>Total footprint area</td><td>{format_number(energy_totals.get("footprint_area_m2", 0.0), 2)} m²</td></tr>
            <tr><td>Total floor area</td><td>{format_number(energy_totals.get("total_floor_area_m2", 0.0), 2)} m²</td></tr>
            <tr><td>Average NPV/m²</td><td>{format_number(energy_totals.get("npv_per_m2_average", 0.0), 2)} SEK/m²</td></tr>
        </table>

        <h3>System CAPEX</h3>
        <table>
            <tr><td>Grid CAPEX</td><td>{format_number(grid_capex_total, 2)} SEK</td></tr>
            <tr><td>DHN CAPEX</td><td>{format_number(dhn_capex_total, 2)} SEK</td></tr>
            <tr><td>Total network CAPEX</td><td>{format_number(grid_capex_total + dhn_capex_total, 2)} SEK</td></tr>
            <tr><td>Total CAPEX incl. buildings</td><td>{format_number(grid_capex_total + dhn_capex_total + energy_totals.get("building_capital_cost", 0.0), 2)} SEK</td></tr>
        </table>

        {change_html}
    </div>

    <style>
        #market-panel {{
            position: fixed;
            top: 10px;
            right: 10px;
            width: 490px;
            max-height: 92vh;
            overflow-y: auto;
            z-index: 9999;
            background: white;
            border: 2px solid #333;
            border-radius: 8px;
            padding: 12px;
            font-family: Arial, sans-serif;
            font-size: 12px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.35);
        }}

        #market-panel h2 {{
            margin-top: 0;
            font-size: 18px;
        }}

        #market-panel h3 {{
            margin-bottom: 5px;
            margin-top: 14px;
            font-size: 14px;
            border-bottom: 1px solid #ccc;
        }}

        #market-panel table {{
            width: 100%;
            border-collapse: collapse;
        }}

        #market-panel td {{
            border-bottom: 1px solid #eee;
            padding: 3px;
            vertical-align: top;
        }}

        #market-panel td:first-child {{
            font-weight: bold;
            width: 62%;
        }}
    </style>
    """

    return html


def npv_color(value, vmin, vmax):
    try:
        value = float(value)
    except Exception:
        return "#cccccc"

    if not np.isfinite(value):
        return "#cccccc"

    if vmax == vmin:
        return "#f2f2f2"

    if value < 0:
        ratio = value / vmin if vmin < 0 else 0.0
        ratio = max(0.0, min(1.0, ratio))

        if ratio > 0.75:
            return "#67001f"
        if ratio > 0.50:
            return "#b2182b"
        if ratio > 0.25:
            return "#ef8a62"
        return "#fddbc7"

    if value > 0:
        ratio = value / vmax if vmax > 0 else 0.0
        ratio = max(0.0, min(1.0, ratio))

        if ratio > 0.75:
            return "#006837"
        if ratio > 0.50:
            return "#1a9850"
        if ratio > 0.25:
            return "#66bd63"
        return "#d9f0d3"

    return "#f7f7f7"


def add_npv_legend(m, vmin, vmax):
    legend_html = f"""
    <div id="npv-legend">
        <h4>NPV / m²</h4>
        <div><span style="background:#67001f"></span> Very negative</div>
        <div><span style="background:#b2182b"></span> Negative</div>
        <div><span style="background:#ef8a62"></span> Slightly negative</div>
        <div><span style="background:#f7f7f7"></span> Around zero</div>
        <div><span style="background:#66bd63"></span> Slightly positive</div>
        <div><span style="background:#1a9850"></span> Positive</div>
        <div><span style="background:#006837"></span> Very positive</div>
        <p>Min: {format_number(vmin, 2)} SEK/m²<br>Max: {format_number(vmax, 2)} SEK/m²</p>
    </div>

    <style>
        #npv-legend {{
            position: fixed;
            left: 10px;
            bottom: 30px;
            z-index: 9999;
            background: white;
            border: 2px solid #333;
            border-radius: 8px;
            padding: 10px;
            font-family: Arial, sans-serif;
            font-size: 12px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.25);
        }}

        #npv-legend h4 {{
            margin: 0 0 6px 0;
            font-size: 13px;
        }}

        #npv-legend span {{
            display: inline-block;
            width: 18px;
            height: 10px;
            margin-right: 6px;
            border: 1px solid #999;
        }}

        #npv-legend p {{
            margin: 6px 0 0 0;
        }}
    </style>
    """

    m.get_root().html.add_child(folium.Element(legend_html))


def add_buildings_to_map(m, buildings_gdf, config=None):
    gdf = as_wgs84(buildings_gdf)

    if gdf is None:
        return

    if "npv_per_m2" in gdf.columns:
        values = pd.to_numeric(gdf["npv_per_m2"], errors="coerce").fillna(0.0)
        vmin = float(values.min())
        vmax = float(values.max())
    else:
        vmin = 0.0
        vmax = 0.0

    def style_function(feature):
        props = feature.get("properties", {})
        value = props.get("npv_per_m2", 0.0)

        return {
            "fillColor": npv_color(value, vmin, vmax),
            "color": "#222222",
            "weight": 1,
            "fillOpacity": 0.70,
        }

    tooltip_fields = [
        c for c in [
            "id",
            "Name",
            "EEdepth",
            "SHSource",
            "DHWsource",
            "PVType",
            "PVpercentage",
            "Floors",
            "footprint_area_m2",
            "total_floor_area_m2",
            "electricity_bought_kwh",
            "electricity_sold_kwh",
            "dhn_bought_kwh",
            "fuel_bought_kwh",
            "fuel_type",
            "primary_energy_kwh",
            "electricity_operational_cost",
            "dhn_operational_cost",
            "fuel_operational_cost",
            "total_operational_cost",
            "building_capital_cost",
            "building_npv",
            "npv_per_m2",
        ]
        if c in gdf.columns
    ]

    fg = folium.FeatureGroup(name="Buildings colored by NPV/m²", show=True)

    folium.GeoJson(
        gdf,
        name="Buildings colored by NPV/m²",
        style_function=style_function,
        tooltip=folium.GeoJsonTooltip(
            fields=tooltip_fields,
            aliases=tooltip_fields,
            localize=True,
        ) if tooltip_fields else None,
    ).add_to(fg)

    fg.add_to(m)
    add_npv_legend(m, vmin, vmax)


def add_lines_to_map(m, line_gdf, name, color, weight=6, show=True):
    gdf = as_wgs84(line_gdf)

    if gdf is None:
        return

    fg = folium.FeatureGroup(name=name, show=show)

    tooltip_fields = [
        c for c in [
            "network_type",
            "dhn_id",
            "grid_id",
            "line_id",
            "start_node",
            "end_node",
            "from_bus",
            "to_bus",
        ]
        if c in gdf.columns
    ]

    folium.GeoJson(
        gdf,
        name=name,
        style_function=lambda feature: {
            "color": color,
            "weight": weight,
            "opacity": 1.0,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=tooltip_fields,
            aliases=tooltip_fields,
        ) if tooltip_fields else None,
    ).add_to(fg)

    fg.add_to(m)


def add_nodes_to_map(m, nodes_gdf):
    gdf = as_wgs84(nodes_gdf)

    if gdf is None:
        return

    dhn_fg = folium.FeatureGroup(name="DHN nodes", show=False)
    grid_fg = folium.FeatureGroup(name="Grid nodes", show=False)

    for _, row in gdf.iterrows():
        y = row.geometry.y
        x = row.geometry.x

        network_type = row.get("network_type", "")
        node_id = row.get("node_id", "")
        node_type = row.get("node_type", "")
        network_id = row.get("network_id", "")

        popup = (
            f"{network_type}<br>"
            f"network_id: {network_id}<br>"
            f"node_id: {node_id}<br>"
            f"node_type: {node_type}"
        )

        if network_type == "DHN node":
            color = "#a50026" if node_type == "supply" else "#d7191c"

            folium.CircleMarker(
                location=[y, x],
                radius=4 if node_type == "supply" else 3,
                color=color,
                fill=True,
                fill_opacity=0.9,
                popup=popup,
            ).add_to(dhn_fg)
        else:
            folium.CircleMarker(
                location=[y, x],
                radius=3,
                color="#005a9c",
                fill=True,
                fill_opacity=0.9,
                popup=popup,
            ).add_to(grid_fg)

    dhn_fg.add_to(m)
    grid_fg.add_to(m)


def fit_map_to_layers(m, *layers):
    all_layers = []

    for layer in layers:
        wgs = as_wgs84(layer)

        if wgs is not None and len(wgs) > 0:
            all_layers.append(wgs)

    if not all_layers:
        return

    merged = pd.concat(all_layers, ignore_index=True)
    merged = gpd.GeoDataFrame(merged, geometry="geometry", crs="EPSG:4326")

    minx, miny, maxx, maxy = merged.total_bounds

    if not np.isfinite([minx, miny, maxx, maxy]).all():
        return

    if minx == maxx or miny == maxy:
        return

    m.fit_bounds([[miny, minx], [maxy, maxx]])


def print_debug_examples(retrofit):
    try:
        dhn_systems = getattr(retrofit, "District_Heating_Systems", None)

        if isinstance(dhn_systems, dict):
            dhn = next(iter(dhn_systems.values()))
        elif dhn_systems:
            dhn = dhn_systems[0]
        else:
            dhn = None

        if dhn is not None:
            if getattr(dhn, "nodes", None):
                print("  DHN node example:", dhn.nodes[0].__dict__)
            if getattr(dhn, "lines", None):
                print("  DHN line example:", dhn.lines[0].__dict__)
            print("  DHN has hourly_heat_demand:", hasattr(dhn, "hourly_heat_demand"))
            print("  DHN has hourly_heat_loss:", hasattr(dhn, "hourly_heat_loss"))
    except Exception as e:
        print("  Could not print DHN examples:", e)

    try:
        grids = getattr(retrofit, "Electrical_Network", None)

        if isinstance(grids, dict):
            grid = next(iter(grids.values()))
        elif grids:
            grid = grids[0]
        else:
            grid = None

        if grid is not None:
            if getattr(grid, "nodes", None):
                print("  Grid node example:", grid.nodes[0].__dict__)
            if getattr(grid, "lines", None):
                print("  Grid line example:", grid.lines[0].__dict__)
    except Exception as e:
        print("  Could not print grid examples:", e)


def write_map(
    game_step,
    buildings_gdf,
    retrofit,
    config,
    optimal_set,
    grid_pricing,
    dhn_pricing,
    capital_costs_grid,
    capital_costs_dhn,
    changes,
    assumptions,
    output_folder,
    open_browser=False,
):
    map_folder = Path(output_folder) / "maps"
    map_folder.mkdir(parents=True, exist_ok=True)

    buildings_map_gdf = attach_energy_attributes_to_buildings(
        buildings_gdf=buildings_gdf,
        optimal_set=optimal_set,
        assumptions=assumptions,
    )

    energy_totals = summarize_map_totals(buildings_map_gdf)

    dhn_gdf = dhn_to_gdf(retrofit, buildings_map_gdf)
    grid_gdf = grid_to_gdf(retrofit, buildings_map_gdf)
    nodes_gdf = nodes_to_gdf(retrofit, buildings_map_gdf)

    dhn_systems = getattr(retrofit, "District_Heating_Systems", None)
    grids = getattr(retrofit, "Electrical_Network", None)

    dhn_operation_totals = compute_dhn_operation_from_objects(
        dhns=dhn_systems,
        unit=assumptions.get("operation_unit", "Wh"),
    )

    source_counts = count_building_heat_sources(
        config=config,
        buildings_gdf=buildings_map_gdf,
    )

    building_count = len(buildings_map_gdf)
    dhn_count = count_container(dhn_systems)
    grid_count = count_container(grids)
    dhn_line_count = 0 if dhn_gdf is None else len(dhn_gdf)
    grid_line_count = 0 if grid_gdf is None else len(grid_gdf)

    print("Map debug:")
    print("  buildings:", building_count, buildings_map_gdf.crs)
    print("  buildings connected to DHN:", source_counts.get("dhn_connected_buildings", 0))
    print("  buildings with heat pump:", source_counts.get("hp_buildings", 0))
    print("  dhn systems:", dhn_count)
    print("  dhn pipes:", dhn_line_count, None if dhn_gdf is None else dhn_gdf.crs)
    print("  grid systems:", grid_count)
    print("  grid cable lines:", grid_line_count, None if grid_gdf is None else grid_gdf.crs)
    print("  dhn operation totals:", dhn_operation_totals)
    print("  energy totals:", energy_totals)

    if dhn_line_count == 0 or grid_line_count == 0:
        print_debug_examples(retrofit)

    center = map_center_from_layers(buildings_map_gdf, dhn_gdf, grid_gdf)

    m = folium.Map(
        location=center,
        zoom_start=16,
        tiles="OpenStreetMap",
        control_scale=True,
    )

    folium.TileLayer("CartoDB positron", name="Light map").add_to(m)
    folium.TileLayer("CartoDB dark_matter", name="Dark map").add_to(m)

    add_buildings_to_map(m, buildings_map_gdf, config=config)
    add_lines_to_map(m, dhn_gdf, name="DHN pipes", color="#d7191c", weight=7, show=True)
    add_lines_to_map(m, grid_gdf, name="Electrical cable lines", color="#2c7bb6", weight=6, show=True)
    add_nodes_to_map(m, nodes_gdf)

    fit_map_to_layers(m, buildings_map_gdf, dhn_gdf, grid_gdf)

    m.add_child(MeasureControl())
    Fullscreen().add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)

    panel_html = make_side_panel_html(
        game_step=game_step,
        capital_costs_grid=capital_costs_grid,
        capital_costs_dhn=capital_costs_dhn,
        changes=changes,
        dhn_count=dhn_count,
        dhn_line_count=dhn_line_count,
        grid_count=grid_count,
        grid_line_count=grid_line_count,
        building_count=building_count,
        energy_totals=energy_totals,
        dhn_operation_totals=dhn_operation_totals,
        source_counts=source_counts,
    )

    m.get_root().html.add_child(folium.Element(panel_html))

    if open_browser:
        m.get_root().header.add_child(
            folium.Element('<meta http-equiv="refresh" content="10">')
        )

    step_path = map_folder / f"game_step_{game_step:03d}.html"
    current_path = map_folder / "current_map.html"

    m.save(step_path)
    shutil.copyfile(step_path, current_path)

    print(f"Map saved: {current_path.resolve()}")

    if open_browser:
        try:
            webbrowser.open(current_path.resolve().as_uri())
        except Exception as e:
            print(f"Could not open browser automatically: {e}")
            print(f"Open manually: {current_path.resolve()}")

    return step_path, current_path