"""
Tests
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os

import pytest
import numpy as np
#########################################################
# Config loading
# Loads a global config object
from eureca_building.config import load_config

config_path = os.path.join('.', 'eureca_building', 'test', 'config.json')
load_config(config_path)
from eureca_building.config import CONFIG

from eureca_building.material import Material, AirGapMaterial
from eureca_building.window import SimpleWindow
from eureca_building.construction import Construction
from eureca_building.construction_dataset import ConstructionDataset
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


class TestMaterials:
    """
    This is a test class for the pytest module.
    It tests Material class and its property
    """

    def test_material(self):
        # Standard material creation
        Material("Test material")

    def test_material_comp(self):
        # Standard material creation
        Material("Test material", thick=0.100, cond=1.00, spec_heat=1000.0, dens=1000.0)

    def test_material_with_prop_wrong(self):
        # Standard material creation
        with pytest.raises(MaterialPropertyOutsideBoundaries):
            Material("Test material", cond=1000.0)

    def test_material_setter(self):
        # Standard material creation
        mat = Material(
            "Test material", thick=0.100, cond=1.00, spec_heat=1000.0, dens=1000.0
        )

        with pytest.raises(MaterialPropertyOutsideBoundaries):
            mat.thick = 100.0

    def test_material_setter_good(self):
        # Standard material creation
        mat = Material(
            "Test material", thick=0.100, cond=1.00, spec_heat=1000.0, dens=1000.0
        )

        mat.dens = 800.0

    def test_material_setter_list(self):
        # Standard material creation
        mat = Material(
            "Test material", thick=0.100, cond=1.00, spec_heat=1000.0, dens=1000.0
        )

        with pytest.raises(TypeError):
            mat.thick = "fd"

    def test_air_material(self):
        AirGapMaterial("Test material")

    def test_air_material_2(self):
        AirGapMaterial("Test material", thick=0.100, thermal_resistance=1)

    def test_airgapmaterial_with_prop_wrong(self):
        with pytest.raises(MaterialPropertyOutsideBoundaries):
            AirGapMaterial("Test material", thermal_resistance=100)

    def test_airgapmaterial_setter(self):
        mat = AirGapMaterial("Test material", thick=0.100, thermal_resistance=1)

        with pytest.raises(MaterialPropertyOutsideBoundaries):
            mat.thick = 100.0

    def test_airgapmaterial_setter_good(self):
        mat = AirGapMaterial("Test material", thick=0.100, thermal_resistance=1)

        mat.thermal_resistance = 2.0

    def test_airgapmaterial_setter_list(self):
        mat = AirGapMaterial("Test material", thick=0.100, thermal_resistance=1)

        with pytest.raises(TypeError):
            mat.thick = "fd"


class TestConstruction:
    """
    This is a test class for the pytest module.
    It tests Construction class and its property
    """

    def test_wall(self):
        plaster = Material(
            "plaster", thick=0.01, cond=1.0, spec_heat=800.0, dens=2000.0
        )

        hollowed_bricks = Material(
            "hollowed_bricks", thick=0.150, cond=1.4, spec_heat=800.0, dens=2000.0,
        )

        air = AirGapMaterial("AirMaterial", thick=0.02, thermal_resistance=0.5)

        insulation = Material("tyles", thick=0.01, cond=1, spec_heat=840.0, dens=2300.0)

        Construction(
            "ExtWall",
            materials_list=[plaster, hollowed_bricks, air, insulation, plaster],
            construction_type="ExtWall",
        )

    def test_wall_values(self):
        plaster = Material(
            "plaster", thick=0.01, cond=1.0, spec_heat=800.0, dens=2000.0
        )

        hollowed_bricks = Material(
            "hollowed_bricks", thick=0.150, cond=1.4, spec_heat=800.0, dens=2000.0,
        )

        air = AirGapMaterial("AirMaterial", thick=0.02, thermal_resistance=0.5)

        insulation = Material(
            "tyles", thick=0.01, cond=0.03, spec_heat=1000.0, dens=30.0
        )

        ext_wall = Construction(
            "ExtWall",
            materials_list=[plaster, hollowed_bricks, air, insulation, plaster],
            construction_type="ExtWall",
        )

        # U net should be 1.041150223
        # U should be 0.884684616

        assert abs(ext_wall._u_value - 0.884684616) < ext_wall._u_value / 0.001
        assert abs(ext_wall._u_value_net - 1.041150223) < ext_wall._u_value_net / 0.001

    def test_construction_from_U_value(self):

        ext_wall = Construction.from_U_value(
            "ExtWall",
            0.7,
            weight_class = "Medium",
            construction_type="ExtWall",
        )

        # U net should be 1.041150223
        # U should be 0.884684616

        assert abs(ext_wall._u_value - 0.7) < ext_wall._u_value / 0.001

    def test_window_values(self):
        window = SimpleWindow(
            name="window_1",
            u_value=5,
            solar_heat_gain_coef=0.2,
            visible_transmittance=0.3,
            frame_factor=0.1,
            shading_coef_int=0.1,
            shading_coef_ext=0.1,
        )


class TestConstructionDataset:
    """
    This is a test class for the pytest module.
    It tests ConstructionDataset class and its property
    """

    def test_read_excel_method(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "materials_and_construction_test.xlsx",
        )

        dataset = ConstructionDataset.read_excel(path)


class TestSurface:
    """
    This is a test class for the pytest module.
    It tests Surface class and its property
    """

    def test_creation_of_surface_zero(self):
        with pytest.raises(SurfaceWrongNumberOfVertices):
            surf = Surface("Surface 1")

    def test_creation_of_surface(self):
        surf = Surface("Surface 1", vertices=((0, 0, 0), (0, 1, 1), (0, 1, 2),))
        print(surf._vertices)

    def test_creation_of_surface_fake(self):
        with pytest.raises(SurfaceWrongNumberOfVertices):
            Surface("Surface 1", vertices=((0, 1, 1), (0, 1, 2),))

    def test_nocoplanar_surface(self):
        with pytest.raises(NonPlanarSurface):
            Surface("Surface 1", vertices=((0, 0, 0), (0, 1, 1), (0, 1, 2), (1, 2, 4)))

    def test_wwr_surface(self):
        surf_1 = Surface(
            "Surface 1", vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)), wwr=0.4,
        )

        assert surf_1._opaque_area == 0.6
        assert surf_1._glazed_area == 0.4

    def test_wwr_surface_2(self):
        with pytest.raises(WindowToWallRatioOutsideBoundaries):
            Surface(
                "Surface 1",
                vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
                wwr=2.0,
            )

    def test_wwr_surface_3(self):
        surf_1 = Surface(
            "Surface 1", vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)), wwr=0.4,
        )

        surf_1._wwr = 0.2

        assert surf_1._opaque_area == 0.8
        assert surf_1._glazed_area == 0.2

    def test_subdivision_solar_calc(self):
        surf_1 = Surface(
            "Surface 1",
            vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
            wwr=0.4,
            subdivisions_solar_calc={
                "azimuth_subdivisions": 3,
                "height_subdivisions": 5,
            },
        )

    def test_subdivision_solar_calc_2(self):
        with pytest.raises(ValueError):
            surf_1 = Surface(
                "Surface 1",
                vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
                wwr=0.4,
                subdivisions_solar_calc={
                    "azimuth_subdivisions": 1000,
                    "height_subdivisions": 5,
                },
            )

    def test_subdivision_solar_calc_3(self):
        with pytest.raises(ValueError):
            surf_1 = Surface(
                "Surface 1",
                vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
                wwr=0.4,
                subdivisions_solar_calc={
                    "azimuth_subdivisions": 3,
                    "height_subdivisions": 500000,
                },
            )

    def test_subdivision_solar_calc_4(self):
        with pytest.raises(TypeError):
            surf_1 = Surface(
                "Surface 1",
                vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
                wwr=0.4,
                subdivisions_solar_calc={
                    "azimuth_subdivisions": "a",
                    "height_subdivisions": 500000,
                },
            )

    def test_subdivision_solar_calc_5(self):
        surf_1 = Surface(
            "Surface 1", vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)), wwr=0.4,
        )
        surf_1.subdivisions_solar_calc = {
            "azimuth_subdivisions": 8,
            "height_subdivisions": 3,
        }

    def test_surface_type(self):
        for s_type in ["ExtWall", "GroundFloor", "Roof"]:
            Surface(
                "Surface 1",
                vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
                wwr=0.4,
                surface_type=s_type,
            )

    def test_surface_type_2(self):
        surf_1 = Surface(
            "Surface 1",
            vertices=((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
            wwr=0.4,
            surface_type="ExtWall",
        )
        with pytest.raises(InvalidSurfaceType):
            surf_1.surface_type = "bds"

    def test_creation_of_surfaceIM_zero(self):
        surf = SurfaceInternalMass("Surface 1")
        print(surf.name)

    def test_creation_of_surfaceIM(self):
        surf = SurfaceInternalMass("Surface 1", area=10.)
        print(surf._area)

    def test_creation_of_surface_fakeIM(self):
        with pytest.raises(NegativeSurfaceArea):
            SurfaceInternalMass("Surface 1", area=-2.)

    def test_surfaceIM_type(self):
        for s_type in ["IntWall", "IntCeiling", "IntFloor"]:
            SurfaceInternalMass(
                "Surface 1",
                surface_type=s_type,
            )

    def test_surfaceIM_type_wrong(self):
        for s_type in ["IntWall1", "IntC", "IFloor"]:
            with pytest.raises(InvalidSurfaceType):
                SurfaceInternalMass(
                    "Surface 1",
                    surface_type=s_type,
                )
