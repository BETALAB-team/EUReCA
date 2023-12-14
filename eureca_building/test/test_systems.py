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

config_path = os.path.join('.', 'eureca_building', 'test', 'config.json')
load_config(config_path)
from eureca_building.config import CONFIG
import pytest
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from eureca_building.material import Material, AirGapMaterial
from eureca_building.window import SimpleWindow
from eureca_building.construction import Construction
from eureca_building.weather import WeatherFile
from eureca_building.systems import hvac_heating_systems_classes, hvac_cooling_systems_classes
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


def plot_contour(x,y,z,resolution = 50,contour_method='linear'):
    resolution = str(resolution)+'j'
    X,Y = np.mgrid[min(x):max(x):complex(resolution),   min(y):max(y):complex(resolution)]
    points = [[a,b] for a,b in zip(x,y)]
    Z = griddata(points, z, (X, Y), method=contour_method)
    return X,Y,Z

class TestSystem:
    """
    This is a test class for the pytest module.
    It tests WeatherFile class and its property
    """

    def test_traditional_boiler(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        weather = WeatherFile(path, time_steps=CONFIG.ts_per_hour)

        heating_system = hvac_heating_systems_classes["TraditionalBoiler"]()
        for pw in range(1,500000,10000):
            heating_system.set_system_capacity(design_power=float(pw), weather = weather)
            for heat in np.linspace(pw*0.2, pw, 5):
                heating_system.solve_system(heat, weather, 50, 20, 50)


    def test_condensing_boiler(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        weather = WeatherFile(path, time_steps=CONFIG.ts_per_hour)

        heating_system = hvac_heating_systems_classes["CondensingBoiler"]()
        for pw in range(1,500000,10000):
            heating_system.set_system_capacity(design_power=float(pw), weather = weather)
            for heat in np.linspace(pw*0.2, pw, 5):
                heating_system.solve_system(heat, weather, 50, 20, 50)

    def test_air_cooler(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        weather = WeatherFile(path, time_steps=CONFIG.ts_per_hour)

        cooling_system = hvac_cooling_systems_classes["SplitAirCooler"]()
        for pw in range(-1,-500000,-10000):
            cooling_system.set_system_capacity(design_power=float(pw), weather = weather)
            for cool in np.linspace(pw*0.2, pw, 5):
                cooling_system.solve_system(cool, weather, 10000, 20, 50)

    def test_air_conditioner(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        weather = WeatherFile(path, time_steps=CONFIG.ts_per_hour)

        cooling_system = hvac_cooling_systems_classes["SplitAirConditioner"]()
        for pw in range(-1,-500000,-10000):
            cooling_system.set_system_capacity(design_power=float(pw), weather = weather)
            for cool in np.linspace(pw*0.2, pw, 5):
                cooling_system.solve_system(cool, weather, 10000, 20, 50)

    def test_chiller_aw(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        weather = WeatherFile(path, time_steps=CONFIG.ts_per_hour)

        cooling_system = hvac_cooling_systems_classes["ChillerAirtoWater"]()
        for pw in range(-1,-500000,-10000):
            cooling_system.set_system_capacity(design_power=float(pw), weather = weather)
            for cool in np.linspace(pw*0.2, pw, 5):
                cooling_system.solve_system(cool, weather, 10000, 20, 50)

    def test_EN15316Sytems(self):
        path = os.path.join(
            "eureca_building",
            "example_scripts",
            "ITA_Venezia-Tessera.161050_IGDG.epw",
        )
        weather = WeatherFile(path, time_steps=CONFIG.ts_per_hour)

        relevant_hvac = [x for x in hvac_heating_systems_classes.keys() if x not in [
            "IdealLoad",
            "CondensingBoiler",
            "TraditionalBoiler"
        ]
                         ]
        for hvac in relevant_hvac:
            heating_system = hvac_heating_systems_classes[hvac](heating_system_key = hvac)
            results = np.array([[0.,0.,0]])
            for pw in range(1,500000,10000):
                heating_system.set_system_capacity(design_power=float(pw), weather = weather)
                for heat in np.linspace(pw*0.2, pw, 5):
                    heating_system.solve_system(heat, weather, 10000, 20, 50)
                    results = np.vstack([results,np.array([pw,heat,heating_system.total_efficiency])])

            # results = results[1:]
            # X, Y, Z = plot_contour(results[:,0], results[:,1], results[:,2], resolution=50, contour_method='linear')
            # fig, [ax1, ax2] = plt.subplots(nrows = 2)
            # # fig = plt.figure(figsize=(13, 6))
            # #
            # # ax1 = fig.add_subplot(121, projection="3d")
            # # ax2 = fig.add_subplot(122, projection="3d")
            # # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.01)
            #
            # if "A-W" in hvac.split(",")[0]:
            #     vmax = 4
            #     vmin = 0
            # else:
            #     vmax = 1.2
            #     vmin = 0.5
            # with plt.style.context("seaborn"):
            #     fig.suptitle(hvac)
            #     img = ax1.contourf(X/1000, Y/1000, Z, cmap='coolwarm',vmin=vmin,vmax=vmax)
            #     fig.colorbar(img)
            #     img = ax2.contourf(X/1000, Y/1000, Z,cmap='coolwarm')
            #     fig.colorbar(img)
            #
            # ax1.set_xlabel("Design load [kW]")
            # ax1.set_xlabel("Partial load [kW]")
            # #ax2.view_init(elev=28, azim=120)
            # # ax2.view_init(elev=5, azim=115)
            # fig.savefig(os.path.join(".","eureca_building","test","hvac_plots",f"{hvac}.png"))
            # plt.close()




        relevant_hvac = [x for x in hvac_cooling_systems_classes.keys() if x not in [
            "IdealLoad",
            "SplitAirCooler",
            "ChillerAirtoWater",
            "SplitAirConditioner"
        ]
                         ]
        for hvac in relevant_hvac:
            results = np.array([[0., 0., 0]])
            cooling_system = hvac_cooling_systems_classes[hvac](cooling_system_key=hvac)
            for pw in range(-1, -500000, -10000):
                cooling_system.set_system_capacity(design_power=float(pw), weather=weather)
                for heat in np.linspace(pw * 0.2, pw, 5):
                    cooling_system.solve_system(heat, weather, 10000, 20, 50)
                    results = np.vstack([results, np.array([pw, heat, cooling_system.total_efficiency])])

            # results = results[1:]
            # X, Y, Z = plot_contour(results[:,0], results[:,1], results[:,2], resolution=50, contour_method='linear')
            # fig, [ax1, ax2] = plt.subplots(nrows = 2)
            # # fig = plt.figure(figsize=(13, 6))
            # #
            # # ax1 = fig.add_subplot(121, projection="3d")
            # # ax2 = fig.add_subplot(122, projection="3d")
            # # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.01)
            #
            # if "A-W" in hvac.split(",")[0]:
            #     vmax = 4
            #     vmin = 0
            # else:
            #     vmax = 1.2
            #     vmin = 0.5
            # with plt.style.context("seaborn"):
            #     fig.suptitle(hvac)
            #     img = ax1.contourf(X/1000, Y/1000, Z, cmap='coolwarm',vmin=vmin,vmax=vmax)
            #     fig.colorbar(img)
            #     img = ax2.contourf(X/1000, Y/1000, Z,cmap='coolwarm')
            #     fig.colorbar(img)
            #
            # ax1.set_xlabel("Design load [kW]")
            # ax1.set_xlabel("Partial load [kW]")
            # #ax2.view_init(elev=28, azim=120)
            # # ax2.view_init(elev=5, azim=115)
            # fig.savefig(os.path.join(".","eureca_building","test","hvac_plots",f"{hvac}.png"))
            # plt.close()
