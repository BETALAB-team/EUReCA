import os
import json
AREA_MAP = {
    "urban": 1.3,      
    "suburban": 0.9,   
    "park": 0.7,       
    "new_dev": 0.5     
}


def build_cost_table(pipe_data):
    """
    Accepts:
    - dict DN->cost
    - JSON-like dict with 'pipes'
    - path to JSON file

    Returns:
    - dict: { "DNxx": cost_per_meter }
    """

    # -------------------------
    # Case 1: string path
    # -------------------------
    if isinstance(pipe_data, str):
        if not os.path.exists(pipe_data):
            raise FileNotFoundError(f"File not found: {pipe_data}")

        with open(pipe_data, "r") as f:
            pipe_data = json.load(f)

    # -------------------------
    # Case 2: already DN dict
    # -------------------------
    if isinstance(pipe_data, dict) and all(
        isinstance(k, str) and k.startswith("DN") for k in pipe_data.keys()
    ):
        return pipe_data

    # -------------------------
    # Case 3: JSON structure
    # -------------------------
    if isinstance(pipe_data, dict) and "pipes" in pipe_data:
        cost_table = {}

        for p in pipe_data["pipes"]:
            dn = p["dn"]
            cost = p.get("cost_per_meter")

            if cost is not None:
                cost_table[dn] = cost

        return cost_table

    raise ValueError("Invalid pipe_data format")


# -------------------------
# NORMALIZE DN FORMAT
# -------------------------
def normalize_dn(dn):
    if isinstance(dn, (int, float)):
        return f"DN{int(dn)}"
    return dn



def normalize_dn(dn):
    """
    Ensures DN is in 'DNxxx' format
    """
    if isinstance(dn, (int, float)):
        return f"DN{int(dn)}"
    return dn


def compute_dhn_cost(
    dhn_pipe_changes,
    area_type,
    assumptions,
    pipe_json
):
    """
    Compute total DHN cost.

    Parameters
    ----------
    dhn_pipe_changes : dict
        {
            "new_pipes": [Line, ...],
            "pipes_changed": [(old_line, new_line), ...]
        }

    area_type : str
        "urban", "suburban", "park", "new_dev"

    assumptions : dict
        {
            "replacement_factor_ground": float,
            "removal_factor": float
        }

    pipe_json : dict
        JSON containing DN and cost_per_meter

    Returns
    -------
    float
        Total cost (SEK)
    """

    if area_type not in AREA_MAP:
        raise ValueError(f"Invalid area_type: {area_type}")

    area_factor = AREA_MAP[area_type]

    cost_table = build_cost_table(pipe_data=pipe_json)

    rg = assumptions["replacement_factor_ground"]
    rf = assumptions["removal_factor"]

    ground_share = 0.55
    rest_share = 0.45

    total_cost = 0.0


    for line in dhn_pipe_changes.get("new_pipes", []):

        dn = normalize_dn(line.dn)

        if dn not in cost_table:
            raise ValueError(f"Missing cost for {dn}")

        L = line.length

        total_cost += cost_table[dn] * L * area_factor

    for old_line, new_line in dhn_pipe_changes.get("pipes_changed", []):

        dn_new = normalize_dn(new_line.dn)
        dn_old = normalize_dn(old_line.dn)

        if dn_new not in cost_table:
            raise ValueError(f"Missing cost for {dn_new}")
        if dn_old not in cost_table:
            raise ValueError(f"Missing cost for {dn_old}")

        L = new_line.length

        C_new = cost_table[dn_new]
        C_old = cost_table[dn_old]

        Cg_new = C_new * ground_share
        Cr_new = C_new * rest_share
        Cr_old = C_old * rest_share

        cost_per_m = (
            rg * Cg_new +     
            1.0 * Cr_new +    
            rf * Cr_old       
        )

        total_cost += cost_per_m * L * area_factor

    return total_cost