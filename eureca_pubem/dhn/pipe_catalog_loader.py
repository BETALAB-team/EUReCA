import json
import numpy as np
from pathlib import Path

_BASE = Path(__file__).parent


with open(_BASE / "_pipe_diameters.json", "r") as f:
    _dn_data = json.load(f)

DN_CATALOG = _dn_data["pipes"]
DN_INNER_D = np.array([p["inner_diameter"] for p in DN_CATALOG])
DN_OUTER_D = np.array([p["outer_diameter"] for p in DN_CATALOG])


with open(_BASE / "_pipe_materials.json", "r") as f:
    _mat_data = json.load(f)

PIPE_MATERIALS = _mat_data["materials"]


with open(_BASE / "_pipe_insulation.json", "r") as f:
    _ins_data = json.load(f)

INSULATION_CLASSES = _ins_data["classes"]
