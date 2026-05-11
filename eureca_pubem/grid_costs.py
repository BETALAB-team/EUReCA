import os
import json


def _load_json_like(data):
    if isinstance(data, str):
        if not os.path.exists(data):
            raise FileNotFoundError(f"File not found: {data}")

        with open(data, "r") as f:
            return json.load(f)

    return data


def _get(obj, key, default=None):
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def _get_length(line):
    for key in ["length", "length_m", "Length", "L"]:
        value = _get(line, key, None)
        if value is not None:
            return float(value)

    geom = _get(line, "geometry", None)
    if geom is not None:
        return float(geom.length)

    raise ValueError("Could not find line length")


def _normalize_key(value):
    if value is None:
        return None
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def _extract_cable_key(line, cable_key=None):
    if cable_key is not None:
        value = _get(line, cable_key, None)
        if value is not None:
            return _normalize_key(value)

    for key in [
        "cable_name",
        "cable",
        "cable_type",
        "cable_id",
        "type",
        "name",
        "id",
        "section_mm2"
    ]:
        value = _get(line, key, None)

        if isinstance(value, dict):
            for subkey in ["name", "id", "section_mm2"]:
                subvalue = value.get(subkey)
                if subvalue is not None:
                    return _normalize_key(subvalue)

        if value is not None:
            return _normalize_key(value)

    raise ValueError("Could not find cable identifier on line")


def build_grid_cost_table(cable_data, area_type=None):
    cable_data = _load_json_like(cable_data)

    if isinstance(cable_data, dict) and all(
        isinstance(k, (str, int, float)) for k in cable_data.keys()
    ) and "cables" not in cable_data:
        return {_normalize_key(k): float(v) for k, v in cable_data.items()}

    if isinstance(cable_data, dict) and "cables" in cable_data:
        cable_data = cable_data["cables"]

    if not isinstance(cable_data, list):
        raise ValueError("Invalid cable_data format")

    cost_table = {}

    for cable in cable_data:
        if area_type is not None:
            area_costs = cable.get("total_cost_per_m_by_area_type", {})
            if area_type not in area_costs:
                raise ValueError(
                    f"Missing area-specific cost for {cable.get('name')}, area_type={area_type}"
                )
            cost = float(area_costs[area_type])
        else:
            cost = float(cable["cost_per_m"])

        possible_keys = [
            cable.get("name"),
            cable.get("id"),
            cable.get("section_mm2")
        ]

        for key in possible_keys:
            if key is not None:
                cost_table[_normalize_key(key)] = cost

    return cost_table


def compute_grid_cost(
    grid_line_changes,
    area_type,
    assumptions,
    cable_json,
    cable_key=None
):
    if area_type not in ["city", "town", "rural_normal"]:
        raise ValueError(f"Invalid area_type: {area_type}")

    cost_table = build_grid_cost_table(
        cable_data=cable_json,
        area_type=area_type
    )

    rg = assumptions.get("replacement_factor_ground", 1.0)
    rf = assumptions.get("removal_factor", 0.0)

    ground_share = assumptions.get("ground_share", 0.55)
    rest_share = assumptions.get("rest_share", 0.45)

    total_cost = 0.0

    for line in grid_line_changes.get("new_lines", []):
        cable = _extract_cable_key(line, cable_key=cable_key)

        if cable not in cost_table:
            raise ValueError(f"Missing cost for cable: {cable}")

        L = _get_length(line)

        total_cost += cost_table[cable] * L

    for old_line, new_line in grid_line_changes.get("lines_changed", []):
        cable_new = _extract_cable_key(new_line, cable_key=cable_key)
        cable_old = _extract_cable_key(old_line, cable_key=cable_key)

        if cable_new not in cost_table:
            raise ValueError(f"Missing cost for new cable: {cable_new}")
        if cable_old not in cost_table:
            raise ValueError(f"Missing cost for old cable: {cable_old}")

        L = _get_length(new_line)

        C_new = cost_table[cable_new]
        C_old = cost_table[cable_old]

        Cg_new = C_new * ground_share
        Cr_new = C_new * rest_share
        Cr_old = C_old * rest_share

        cost_per_m = (
            rg * Cg_new +
            Cr_new +
            rf * Cr_old
        )

        total_cost += cost_per_m * L

    return total_cost