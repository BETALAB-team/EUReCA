"""
This module includes a ConstructionDataset, a container class for materials, constructions, and windows
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import pandas as pd
import numpy as np

from eureca_building.construction import Construction
from eureca_building.window import SimpleWindow
from eureca_building.material import Material, AirGapMaterial
from eureca_building.exceptions import WrongConstructionType, WrongMaterialType


class ConstructionDataset:
    """This class is a class to include the list of materials,
    construction and windows that are used in the project

    Attributes
    ----------
    materials_dict : dict
        Dict: dictionary with all Materials
    constructions_dict : dict
        Dict: dictionary with all Constructions
    windows_dict : dict
        Dict: dictionary with all Windows
    """

    def __init__(self):
        """Generates the ConstructionDataset and the list of materials,
        contructions, and windows
        """

        self.materials_dict = {}
        self.constructions_dict = {}
        self.windows_dict = {}

    @classmethod
    def read_excel(cls, file):
        """Read the Materials, Windows and Constrcutions from a spreadsheet.
        See and example at https://github.com/BETALAB-team/EUReCA/tree/main/eureca_building/example_scripts under materials_and_constructions.xlsx

        Parameters
        ----------
        file : str
            path to the file.xls

        Returns
        -------
        eureca_building.construction_dataset.ConstructionDataset
            the object with the three dictionaries including all Materials Windows and Constructions from the spreadsheet
        Notes
        -----
        The Excel file must contain the following sheets:
            - "Materials"
            - "Windows"
            - "Constructions"
            
        Each sheet must follow the naming and units expected by the parser.

        """
        dataset = cls()

        data = pd.read_excel(file, sheet_name=None, index_col=0)

        # Windows
        for win_idx in data["Windows"].index:
            win = data["Windows"].loc[win_idx]
            dataset.windows_dict[win_idx] = SimpleWindow(
                name=win["name"],
                u_value=win["U [W/(m²K)]"],
                solar_heat_gain_coef=win["Solar_Heat_Gain_Coef [-]"],
                visible_transmittance=win["visible_transmittance [-]"],
                frame_factor=win["frame_factor [-]"],
                shading_coef_int=win["shading_coef_int [-]"],
                shading_coef_ext=win["shading_coef_ext [-]"],
            )
        # Materials
        for mat_idx in data["Materials"].index:
            mat = data["Materials"].loc[mat_idx]
            if mat["Material_type"] == "Opaque":
                dataset.materials_dict[mat_idx] = Material(
                    name=mat["name"],
                    thick=mat["Thickness [m]"],
                    cond=mat["Conductivity [W/(m K)]"],
                    spec_heat=mat["Specific_heat [J/(kg K)]"],
                    dens=mat["Density [kg/m³]"],
                )
            elif mat["Material_type"] == "AirGapMaterial":
                dataset.materials_dict[mat_idx] = AirGapMaterial(
                    name=mat["name"],
                    thermal_resistance=mat["Thermal_resistance [(m2 K)/W]"],
                )
            else:
                raise WrongMaterialType(
                    f"Material {mat['name']}, invalid material type"
                )
        # Constructions
        for cons_idx in data["Constructions"].index:
            cons = data["Constructions"].loc[cons_idx]
            target_columns = [f'material_{i}_id' for i in range(1, 11)]
            list_of_materials = [
                dataset.materials_dict[x]
                for col in target_columns
                if col in cons and str(cons[col]) != "nan"
                for x in [cons[col]]]
            
            if (
                    isinstance((cons["U-value [W/(m2 K)]"]), (int, float)) and
                    (cons["Weight class"] in ["Very heavy", "Heavy", "Medium", "Light", "Very light"])
            ):
                dataset.constructions_dict[cons_idx] = Construction.from_U_value(
                                                            cons["name"],
                                                            u_value = float(cons["U-value [W/(m2 K)]"]),
                                                            weight_class =  cons["Weight class"],
                                                            construction_type = cons["type"])

            elif(len(list_of_materials)>0):

                dataset.constructions_dict[cons_idx] = Construction(
                    cons["name"],
                    materials_list=list_of_materials,
                    construction_type=cons["type"],
                )
        return dataset
