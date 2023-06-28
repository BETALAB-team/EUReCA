"""
Tests
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os
from eureca_building.config import load_config

config_path = os.path.join('.', 'eureca_building', 'test', 'config.ini')
load_config(config_path)
from eureca_building.config import CONFIG
import pytest
import numpy as np

from eureca_building.material import Material, AirGapMaterial
from eureca_building.window import SimpleWindow
from eureca_building.construction import Construction
from eureca_building.weather import WeatherFile
from eureca_building.construction_dataset import ConstructionDataset
from eureca_building.thermal_zone import ThermalZone
from eureca_building.surface import Surface, SurfaceInternalMass
from eureca_building.exceptions import (
    MaterialPropertyOutsideBoundaries,
    MaterialPropertyNotFound,
    WrongConstructionType,
    Non3ComponentsVertex,
    SurfaceWrongNumberOfVertices,
    WindowToWallRatioOutsideBoundaries,
    InvalidSurfaceType,
    NonPlanarSurface,
    NegativeSurfaceArea,
)


class TestWeatherFile:
    """
    This is a test class for the pytest module.
    It tests WeatherFile class and its property
    """

    def test_weather_file(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        WeatherFile(path, time_steps=4)

    def test_weather_file_2(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        WeatherFile(path, azimuth_subdivisions=6, height_subdivisions=2)


class TestThermalZone:
    """
    This is a test class for the pytest module.
    It tests ThermalZone class and its property
    """

    def test_thermal_zone(self):
        # Load the material and contructions
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "materials_and_construction_test.xlsx",
        )
        # Define some constructions
        dataset = ConstructionDataset.read_excel(path)
        ceiling = dataset.constructions_dict[13]
        floor = dataset.constructions_dict[13]
        ext_wall = dataset.constructions_dict[35]
        window = dataset.windows_dict[2]

        # Definition of surfaces
        wall_1 = Surface(
            "Wall 1",
            vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
            wwr=0.4,
            surface_type="ExtWall",
            construction=ext_wall,
            window=window,
        )
        wall_2 = Surface(
            "Wall 2",
            vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
            surface_type="ExtWall",
            construction=ext_wall,
        )
        floor_1 = Surface(
            "Fllor 1",
            vertices=((0, 0, 0), (0, 1, 0), (1, 1, 0), (0, 1, 0)),
            wwr=0.0,
            surface_type="GroundFloor",
            construction=floor,
        )
        roof_1 = Surface(
            "Roof 1",
            vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
            wwr=0.4,
            surface_type="Roof",
            construction=ceiling,
            window=window,
        )
        intwall_1 = SurfaceInternalMass(
            "Roof 1",
            area=2.,
            surface_type="IntWall",
            construction=ext_wall,
        )

        # Create zone
        tz1 = ThermalZone(
            name="Zone 1",
            surface_list=[wall_1, wall_2, floor_1, roof_1, intwall_1],
            net_floor_area=None,
            volume=None)

        tz1._ISO13790_params()
        tz1._VDI6007_params()
