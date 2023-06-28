"""
This module includes classes and functions to solve the window system
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import pandas as pd
import numpy as np
from scipy.interpolate import splrep

from eureca_building.exceptions import MaterialPropertyOutsideBoundaries
from eureca_building.units import units, window_material_limits


#%% Simple Window class: refer to EnergyPlus reference


class SimpleWindow(object):

    """
    Defines the simple window model with all its characteristics
    and apply the simple glazing model to estimates the curve of SHGC(theta)
    
    __init__: 
        list of properties: list of the properties of the window

    simpleGlazingModel generates the profiles of SHGC. No input
    
    Methods:
        init
        simpleGlazingModel
    """

    def __init__(
        self,
        name: str,
        u_value: float,
        solar_heat_gain_coef: float,
        visible_transmittance: float = 0.9,
        frame_factor: float = 0.9,
        shading_coef_int: float = 0.05,
        shading_coef_ext: float = 0.05,
    ):
        """
        init method. Stores the data

        Parameters
        ----------
        name : str
            name.
        u_value : float
            u_value [W/(m2 K)].
        solar_heat_gain_coef : float
            [-].
        visible_transmittance : float, optional
            [-]. The default is 0.9.
        frame_factor : float, optional
            [-]. The default is 0.9.
        shading_coef_int : float, optional
            [-]. The default is 0.05.
        shading_coef_ext : float, optional
            [-]. The default is 0.05.

        Returns
        -------
        None.

        """

        self.name = name
        self.u_value = u_value
        self.solar_heat_gain_coef = solar_heat_gain_coef
        self.visible_transmittance = visible_transmittance
        self.frame_factor = frame_factor
        self.shading_coef_int = shading_coef_int
        self.shading_coef_ext = shading_coef_ext

        # Runs the simpleGlazingModel
        self.simpleGlazingModel()

    @property
    def u_value(self) -> float:
        return self._u_value

    @u_value.setter
    def u_value(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Material {self.name}, u_value is not a float: {value}")
        if (
            value < window_material_limits["window_u_value"][0]
            or value > window_material_limits["window_u_value"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "U_value",
                lim=window_material_limits["window_u_value"],
                unit=units["U_value"],
                value=value,
            )
        self._u_value = value

    @property
    def solar_heat_gain_coef(self) -> float:
        return self._solar_heat_gain_coef

    @solar_heat_gain_coef.setter
    def solar_heat_gain_coef(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, solar_heat_gain_coef is not a float: {value}"
            )
        if (
            value < window_material_limits["solar_heat_gain_coefficient"][0]
            or value > window_material_limits["solar_heat_gain_coefficient"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "solar_heat_gain_coefficient",
                lim=window_material_limits["solar_heat_gain_coefficient"],
                unit=units["solar_heat_gain_coefficient"],
                value=value,
            )
        self._solar_heat_gain_coef = value

    @property
    def visible_transmittance(self) -> float:
        return self._visible_transmittance

    @visible_transmittance.setter
    def visible_transmittance(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, visible_transmittance is not a float: {value}"
            )
        if (
            value < window_material_limits["non_dimensional_coefficient"][0]
            or value > window_material_limits["non_dimensional_coefficient"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "visible_transiìmittance",
                lim=window_material_limits["non_dimensional_coefficient"],
                unit=units["non_dimensional_coefficient"],
                value=value,
            )
        self._visible_transmittance = value

    @property
    def frame_factor(self) -> float:
        return self._frame_factor

    @frame_factor.setter
    def frame_factor(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, frame_factor is not a float: {value}"
            )
        if (
            value < window_material_limits["non_dimensional_coefficient"][0]
            or value > window_material_limits["non_dimensional_coefficient"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "frame_factor",
                lim=window_material_limits["non_dimensional_coefficient"],
                unit=units["non_dimensional_coefficient"],
                value=value,
            )
        self._frame_factor = value

    @property
    def shading_coef_int(self) -> float:
        return self._shading_coef_int

    @shading_coef_int.setter
    def shading_coef_int(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, shading_coef_int is not a float: {value}"
            )
        if (
            value < window_material_limits["non_dimensional_coefficient"][0]
            or value > window_material_limits["non_dimensional_coefficient"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "shading_coef_int",
                lim=window_material_limits["non_dimensional_coefficient"],
                unit=units["non_dimensional_coefficient"],
                value=value,
            )
        self._shading_coef_int = value

    @property
    def shading_coef_ext(self) -> float:
        return self._shading_coef_int

    @shading_coef_ext.setter
    def shading_coef_ext(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, shading_coef_ext is not a float: {value}"
            )
        if (
            value < window_material_limits["non_dimensional_coefficient"][0]
            or value > window_material_limits["non_dimensional_coefficient"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "shading_coef_ext",
                lim=window_material_limits["non_dimensional_coefficient"],
                unit=units["non_dimensional_coefficient"],
                value=value,
            )
        self._shading_coef_ext = value

    def simpleGlazingModel(self):

        """
        Description
        This function implement the simplified model for glazing that is used in energy plus
        The output variables are:
        -   4 vectors containing solar transmittance, reflectance, absorptance, SHGC 
            in function of the incident angle  from 0� to 90�, with a 10� step 
        -   the solar transmittance, reflectance and absorption at normal incident;
        -   the equivalent conductivity [W/m K] of the glazing material;
        -   the equivalent thickness [m] of the glazing layer;
        -   the visible reflectance on the front side and back side;
        -   the inward flowing fraction;
        
        ---------ATTENTION-----------
        For some combinations of variable SHGC and U the curves are not
        physically possible,(reflectance > 1). 
        For this reason it is necessary to manipulate the final reflectance and transmittance 
        vectors to manage the curves.
        It is also necessary to force the absolute reflectance to be 1 and the
        absolute transmittance to be 0 for the perpendicular incidence.
        
        
        Controlling input variable
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None        
        """

        if self.u_value < 0:
            print("Negative window U-value. Simulation won" "t proceed.")
            return
        if self.u_value > 7:
            print(
                "U-value may be to high, the model could evaluate un-proper output variable."
            )
        if self.solar_heat_gain_coef < 0 or self.solar_heat_gain_coef > 1:
            print("Solar gain of the window not allowed. Simulation won" "t proceed.")
            return
        # Calculation of thermal resistances W/m^2 K

        if self.u_value < 5.85:
            self.Ri_w = 1 / (0.359073 * np.log(self.u_value) + 6.949915)
        else:
            self.Ri_w = 1 / (1.788041 * self.u_value - 2.886625)
        self.Ro_w = 1 / (0.025342 * self.u_value + 29.163853)
        self.Rl_w = 1 / self.u_value - self.Ro_w - self.Ri_w

        if 1 / self.Rl_w > 7:  # Thickness d of the equivalent layer, m
            self.d = 0.002
        else:
            self.d = 0.05914 - 0.00714 / self.Rl_w
        self.Keff = (
            self.d / self.Rl_w
        )  # Thermal conductivity of the representative layer, W/m K

        # Calculation of solar transmittance for normal incident

        if self.u_value > 4.5:
            if self.solar_heat_gain_coef < 0.7206:
                self.Ts = (0.939998 * self.solar_heat_gain_coef ** 2) + (
                    0.20332 * self.solar_heat_gain_coef
                )
            else:
                self.Ts = (1.30415 * self.solar_heat_gain_coef) - 0.30515
        elif self.u_value < 3.4:
            if self.solar_heat_gain_coef <= 0.15:
                self.Ts = 0.41040 * self.solar_heat_gain_coef
            else:
                self.Ts = (
                    (0.085775 * self.solar_heat_gain_coef ** 2)
                    + (0.963954 * self.solar_heat_gain_coef)
                    - 0.084958
                )
        else:
            if self.solar_heat_gain_coef <= 0.15:
                self.Ts_1 = (0.939998 * self.solar_heat_gain_coef ** 2) + (
                    0.20332 * self.solar_heat_gain_coef
                )
                self.Ts_2 = 0.41040 * self.solar_heat_gain_coef
                self.Ts = self.Ts_2 + (self.Ts_1 - self.Ts_2) / (4.5 - 3.4) * (
                    self.u_value - 3.4
                )
            elif (
                self.solar_heat_gain_coef > 0.15 and self.solar_heat_gain_coef < 0.7206
            ):
                self.Ts_1 = (0.939998 * self.solar_heat_gain_coef ** 2) + (
                    0.20332 * self.solar_heat_gain_coef
                )
                self.Ts_2 = (
                    (0.085775 * self.solar_heat_gain_coef ** 2)
                    + (0.963954 * self.solar_heat_gain_coef)
                    - 0.084958
                )
                self.Ts = self.Ts_2 + (self.Ts_1 - self.Ts_2) / (4.5 - 3.4) * (
                    self.u_value - 3.4
                )
            else:
                self.Ts_1 = (1.30415 * self.solar_heat_gain_coef) - 0.30515
                self.Ts_2 = (
                    (0.085775 * self.solar_heat_gain_coef ** 2)
                    + (0.963954 * self.solar_heat_gain_coef)
                    - 0.084958
                )
                self.Ts = self.Ts_2 + (self.Ts_1 - self.Ts_2) / (4.5 - 3.4) * (
                    self.u_value - 3.4
                )
        # Calculation of inside and outside film resistances in summer condition

        self.x = self.solar_heat_gain_coef - self.Ts

        if self.u_value > 4.5:
            self.Ri_s = 1 / (
                (29.436546 * self.x ** 3)
                - (21.943415 * self.x ** 2)
                + (9.945872 * self.x)
                + 7.426151
            )
            self.Ro_s = 1 / ((2.225824 * self.x) + 20.577080)
        elif self.u_value < 3.4:
            self.Ri_s = 1 / (
                (199.8208128 * self.x ** 3)
                - (90.639733 * self.x ** 2)
                + (19.737055 * self.x)
                + 6.766575
            )
            self.Ro_s = 1 / ((4.475553 * self.x) + 20.674424)
        else:
            self.Ri_s_1 = 1 / (
                (29.436546 * self.x ** 3)
                - (21.943415 * self.x ** 2)
                + (9.945872 * self.x)
                + 7.426151
            )
            self.Ri_s_2 = 1 / (
                (199.8208128 * self.x ** 3)
                - (90.639733 * self.x ** 2)
                + (19.737055 * self.x)
                + 6.766575
            )
            self.Ri_s = self.Ri_s_2 + (self.Ri_s_1 - self.Ri_s_2) / (4.5 - 3.4) * (
                self.u_value - 3.4
            )

            self.Ro_s_1 = 1 / ((2.225824 * self.x) + 20.577080)
            self.Ro_s_2 = 1 / ((4.475553 * self.x) + 20.674424)
            self.Ro_s = self.Ro_s_2 + (self.Ro_s_1 - self.Ro_s_2) / (4.5 - 3.4) * (
                self.u_value - 3.4
            )
        # Inward flowing fraction
        self.Rl = self.Rl_w
        self.N = (self.Ro_s + 0.5 * self.Rl) / (self.Ro_s + self.Rl + self.Ri_s)

        self.As = self.x / self.N
        self.Rs_f = 1 - self.Ts - self.As
        self.Rs_b = self.Rs_f

        # Evaluating the visible proprties of the equivalent layer

        self.Rv_f = (
            -0.0622 * self.visible_transmittance ** 3
            + 0.4277 * self.visible_transmittance ** 2
            - 0.4169 * self.visible_transmittance
            + 0.2399
        )
        self.Rv_b = (
            -0.7409 * self.visible_transmittance ** 3
            + 1.6531 * self.visible_transmittance ** 2
            - 1.2299 * self.visible_transmittance
            + 0.4545
        )

        # Definition of the polinomials
        # Loading polynomial coefficients

        TransCoef = [
            [0.014700000, 1.486000000, -3.852000000, 3.355000000, -0.001474000],
            [0.554600000, 0.035630000, -2.416000000, 2.831000000, -0.002037000],
            [0.770900000, -0.638300000, -1.576000000, 2.448000000, -0.002042000],
            [0.346200000, 0.396300000, -2.582000000, 2.845000000, -0.000280400],
            [2.883000000, -5.873000000, 2.489000000, 1.510000000, -0.002577000],
            [3.025000000, -6.366000000, 3.137000000, 1.213000000, -0.001367000],
            [3.229000000, -6.844000000, 3.535000000, 1.088000000, -0.002891000],
            [3.334000000, -7.131000000, 3.829000000, 0.976600000, -0.002952000],
            [3.146000000, -6.855000000, 3.931000000, 0.786000000, -0.002934000],
            [3.744000000, -8.836000000, 6.018000000, 0.084070000, 0.000482500],
        ]

        RefCoef = [
            [16.320000, -57.820000, 79.240000, -50.080000, 13.340000],
            [40.480000, -119.300000, 134.800000, -70.970000, 16.110000],
            [57.490000, -164.500000, 178.000000, -88.750000, 18.840000],
            [5.714000, -16.670000, 18.630000, -9.756000, 3.074000],
            [-0.548800, -6.498000, 21.200000, -20.970000, 7.814000],
            [4.290000, -12.670000, 14.660000, -8.153000, 2.871000],
            [21.740000, -64.440000, 74.890000, -41.790000, 10.620000],
            [4.341000, -12.800000, 14.780000, -8.203000, 2.879000],
            [41.360000, -117.800000, 127.600000, -64.370000, 14.260000],
            [4.490000, -12.660000, 13.970000, -7.501000, 2.693000],
        ]

        # Defining the 10 correlations
        Ts_A = (
            lambda Cos: TransCoef[0][0] * Cos ** 4
            + TransCoef[0][1] * Cos ** 3
            + TransCoef[0][2] * Cos ** 2
            + TransCoef[0][3] * Cos
            + TransCoef[0][4]
        )
        Ts_B = (
            lambda Cos: TransCoef[1][0] * Cos ** 4
            + TransCoef[1][1] * Cos ** 3
            + TransCoef[1][2] * Cos ** 2
            + TransCoef[1][3] * Cos
            + TransCoef[1][4]
        )
        Ts_C = (
            lambda Cos: TransCoef[2][0] * Cos ** 4
            + TransCoef[2][1] * Cos ** 3
            + TransCoef[2][2] * Cos ** 2
            + TransCoef[2][3] * Cos
            + TransCoef[2][4]
        )
        Ts_D = (
            lambda Cos: TransCoef[3][0] * Cos ** 4
            + TransCoef[3][1] * Cos ** 3
            + TransCoef[3][2] * Cos ** 2
            + TransCoef[3][3] * Cos
            + TransCoef[3][4]
        )
        Ts_E = (
            lambda Cos: TransCoef[4][0] * Cos ** 4
            + TransCoef[4][1] * Cos ** 3
            + TransCoef[4][2] * Cos ** 2
            + TransCoef[4][3] * Cos
            + TransCoef[4][4]
        )
        Ts_F = (
            lambda Cos: TransCoef[5][0] * Cos ** 4
            + TransCoef[5][1] * Cos ** 3
            + TransCoef[5][2] * Cos ** 2
            + TransCoef[5][3] * Cos
            + TransCoef[5][4]
        )
        Ts_G = (
            lambda Cos: TransCoef[6][0] * Cos ** 4
            + TransCoef[6][1] * Cos ** 3
            + TransCoef[6][2] * Cos ** 2
            + TransCoef[6][3] * Cos
            + TransCoef[6][4]
        )
        Ts_H = (
            lambda Cos: TransCoef[7][0] * Cos ** 4
            + TransCoef[7][1] * Cos ** 3
            + TransCoef[7][2] * Cos ** 2
            + TransCoef[7][3] * Cos
            + TransCoef[7][4]
        )
        Ts_I = (
            lambda Cos: TransCoef[8][0] * Cos ** 4
            + TransCoef[8][1] * Cos ** 3
            + TransCoef[8][2] * Cos ** 2
            + TransCoef[8][3] * Cos
            + TransCoef[8][4]
        )
        Ts_J = (
            lambda Cos: TransCoef[9][0] * Cos ** 4
            + TransCoef[9][1] * Cos ** 3
            + TransCoef[9][2] * Cos ** 2
            + TransCoef[9][3] * Cos
            + TransCoef[9][4]
        )

        Rs_A = (
            lambda Cos: RefCoef[0][0] * Cos ** 4
            + RefCoef[0][1] * Cos ** 3
            + RefCoef[0][2] * Cos ** 2
            + RefCoef[0][3] * Cos
            + RefCoef[0][4]
        )
        Rs_B = (
            lambda Cos: RefCoef[1][0] * Cos ** 4
            + RefCoef[1][1] * Cos ** 3
            + RefCoef[1][2] * Cos ** 2
            + RefCoef[1][3] * Cos
            + RefCoef[1][4]
        )
        Rs_C = (
            lambda Cos: RefCoef[2][0] * Cos ** 4
            + RefCoef[2][1] * Cos ** 3
            + RefCoef[2][2] * Cos ** 2
            + RefCoef[2][3] * Cos
            + RefCoef[2][4]
        )
        Rs_D = (
            lambda Cos: RefCoef[3][0] * Cos ** 4
            + RefCoef[3][1] * Cos ** 3
            + RefCoef[3][2] * Cos ** 2
            + RefCoef[3][3] * Cos
            + RefCoef[3][4]
        )
        Rs_E = (
            lambda Cos: RefCoef[4][0] * Cos ** 4
            + RefCoef[4][1] * Cos ** 3
            + RefCoef[4][2] * Cos ** 2
            + RefCoef[4][3] * Cos
            + RefCoef[4][4]
        )
        Rs_F = (
            lambda Cos: RefCoef[5][0] * Cos ** 4
            + RefCoef[5][1] * Cos ** 3
            + RefCoef[5][2] * Cos ** 2
            + RefCoef[5][3] * Cos
            + RefCoef[5][4]
        )
        Rs_G = (
            lambda Cos: RefCoef[6][0] * Cos ** 4
            + RefCoef[6][1] * Cos ** 3
            + RefCoef[6][2] * Cos ** 2
            + RefCoef[6][3] * Cos
            + RefCoef[6][4]
        )
        Rs_H = (
            lambda Cos: RefCoef[7][0] * Cos ** 4
            + RefCoef[7][1] * Cos ** 3
            + RefCoef[7][2] * Cos ** 2
            + RefCoef[7][3] * Cos
            + RefCoef[7][4]
        )
        Rs_I = (
            lambda Cos: RefCoef[8][0] * Cos ** 4
            + RefCoef[8][1] * Cos ** 3
            + RefCoef[8][2] * Cos ** 2
            + RefCoef[8][3] * Cos
            + RefCoef[8][4]
        )
        Rs_J = (
            lambda Cos: RefCoef[9][0] * Cos ** 4
            + RefCoef[9][1] * Cos ** 3
            + RefCoef[9][2] * Cos ** 2
            + RefCoef[9][3] * Cos
            + RefCoef[9][4]
        )

        # Evaluating curves for each zone

        alpha = np.linspace(10, 90, 9)
        Cosalpha = np.cos(np.deg2rad(alpha))

        if self.u_value <= 1.42 and self.solar_heat_gain_coef > 0.45:  # Zone 1
            Ts_alpha = Ts_E(Cosalpha)
            Rs_alpha = Rs_E(Cosalpha)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value <= 1.42
            and self.solar_heat_gain_coef <= 0.45
            and self.solar_heat_gain_coef > 0.35
        ):  # Zone 2 linear interpolation
            Ts1 = Ts_J(Cosalpha)
            Rs1 = Rs_J(Cosalpha)
            Ts2 = Ts_E(Cosalpha)
            Rs2 = Rs_E(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.solar_heat_gain_coef - 0.35) / (
                0.45 - 0.35
            )
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.solar_heat_gain_coef - 0.35) / (
                0.45 - 0.35
            )
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif self.u_value <= 1.42 and self.solar_heat_gain_coef <= 0.35:  # Zone 3
            Ts_alpha = Ts_J(Cosalpha)
            Rs_alpha = Rs_J(Cosalpha)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 1.42
            and self.u_value <= 1.7
            and self.solar_heat_gain_coef > 0.55
        ):  # Zone 4 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = Ts_E(Cosalpha)
            Rs2 = Rs_E(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 1.42
            and self.u_value <= 1.7
            and self.solar_heat_gain_coef <= 0.55
            and self.solar_heat_gain_coef > 0.5
        ):  # Zone 5 bilinear interpolation
            Ts11 = Ts_E(Cosalpha)
            Rs11 = Rs_E(Cosalpha)
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs21 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts22 = Ts_E(Cosalpha)
            Rs22 = Rs_E(Cosalpha)
            Ts_alpha = (
                Ts11 * (0.55 - self.solar_heat_gain_coef) * (1.7 - self.u_value)
                + Ts21 * (self.solar_heat_gain_coef - 0.5) * (1.7 - self.u_value)
                + Ts12 * (0.55 - self.solar_heat_gain_coef) * (self.u_value - 1.42)
                + Ts22 * (self.solar_heat_gain_coef - 0.5) * (self.u_value - 1.42)
            ) / ((0.55 - 0.5) * (1.7 - 1.42))
            Rs_alpha = (
                Rs11 * (0.55 - self.solar_heat_gain_coef) * (1.7 - self.u_value)
                + Rs21 * (self.solar_heat_gain_coef - 0.5) * (1.7 - self.u_value)
                + Rs12 * (0.55 - self.solar_heat_gain_coef) * (self.u_value - 1.42)
                + Rs22 * (self.solar_heat_gain_coef - 0.5) * (self.u_value - 1.42)
            ) / ((0.55 - 0.5) * (1.7 - 1.42))
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 1.42
            and self.u_value <= 1.7
            and self.solar_heat_gain_coef <= 0.5
            and self.solar_heat_gain_coef > 0.45
        ):  # Zone 6 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs2 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 1.42
            and self.u_value <= 1.7
            and self.solar_heat_gain_coef <= 0.45
            and self.solar_heat_gain_coef > 0.35
        ):  # Zone 7 bilinear interpolation
            Ts11 = Ts_J(Cosalpha)
            Rs11 = Rs_J(Cosalpha)
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs21 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts22 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs22 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts_alpha = (
                Ts11 * (0.45 - self.solar_heat_gain_coef) * (1.7 - self.u_value)
                + Ts21 * (self.solar_heat_gain_coef - 0.35) * (1.7 - self.u_value)
                + Ts12 * (0.45 - self.solar_heat_gain_coef) * (self.u_value - 1.42)
                + Ts22 * (self.solar_heat_gain_coef - 0.35) * (self.u_value - 1.42)
            ) / ((0.45 - 0.35) * (1.7 - 1.42))
            Rs_alpha = (
                Rs11 * (0.45 - self.solar_heat_gain_coef) * (1.7 - self.u_value)
                + Rs21 * (self.solar_heat_gain_coef - 0.35) * (1.7 - self.u_value)
                + Rs12 * (0.45 - self.solar_heat_gain_coef) * (self.u_value - 1.42)
                + Rs22 * (self.solar_heat_gain_coef - 0.35) * (self.u_value - 1.42)
            ) / ((0.45 - 0.35) * (1.7 - 1.42))
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 1.42
            and self.u_value <= 1.7
            and self.solar_heat_gain_coef <= 0.35
            and self.solar_heat_gain_coef > 0.3
        ):  # Zone 8 linear interpolation
            Ts1 = Ts_J(Cosalpha)
            Rs1 = Rs_J(Cosalpha)
            Ts2 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs2 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 1.42
            and self.u_value <= 1.7
            and self.solar_heat_gain_coef <= 0.3
            and self.solar_heat_gain_coef > 0.25
        ):  # Zone 9 bilinear interpolation
            Ts11 = Ts_J(Cosalpha)
            Rs11 = Rs_J(Cosalpha)
            Ts12 = Ts_J(Cosalpha)
            Rs12 = Rs_J(Cosalpha)
            Ts21 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)], axis=0)
            Rs21 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)], axis=0)
            Ts22 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs22 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts_alpha = (
                Ts11 * (0.30 - self.solar_heat_gain_coef) * (1.7 - self.u_value)
                + Ts21 * (self.solar_heat_gain_coef - 0.25) * (1.7 - self.u_value)
                + Ts12 * (0.30 - self.solar_heat_gain_coef) * (self.u_value - 1.42)
                + Ts22 * (self.solar_heat_gain_coef - 0.25) * (self.u_value - 1.42)
            ) / ((0.30 - 0.25) * (1.7 - 1.42))
            Rs_alpha = (
                Rs11 * (0.30 - self.solar_heat_gain_coef) * (1.7 - self.u_value)
                + Rs21 * (self.solar_heat_gain_coef - 0.25) * (1.7 - self.u_value)
                + Rs12 * (0.30 - self.solar_heat_gain_coef) * (self.u_value - 1.42)
                + Rs22 * (self.solar_heat_gain_coef - 0.25) * (self.u_value - 1.42)
            ) / ((0.30 - 0.25) * (1.7 - 1.42))
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 1.42
            and self.u_value <= 1.7
            and self.solar_heat_gain_coef <= 0.25
        ):  # Zone 10 linear interpolation
            Ts1 = Ts_J(Cosalpha)
            Rs1 = Rs_J(Cosalpha)
            Ts2 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)], axis=0)
            Rs2 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)], axis=0)
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 1.42) / (1.7 - 1.42)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value <= 3.41
            and self.u_value > 1.7
            and self.solar_heat_gain_coef > 0.55
        ):  # Zone 11
            Ts_alpha = Ts_E(Cosalpha)
            Rs_alpha = Rs_E(Cosalpha)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value <= 3.41
            and self.u_value > 1.7
            and self.solar_heat_gain_coef <= 0.55
            and self.solar_heat_gain_coef > 0.5
        ):  # Zone 12
            Ts1 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs1 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts2 = Ts_E(Cosalpha)
            Rs2 = Rs_E(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.solar_heat_gain_coef - 0.5) / (
                0.55 - 0.5
            )
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.solar_heat_gain_coef - 0.5) / (
                0.55 - 0.5
            )
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value <= 3.41
            and self.u_value > 1.7
            and self.solar_heat_gain_coef <= 0.5
            and self.solar_heat_gain_coef > 0.3
        ):  # Zone 13
            Ts_alpha = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs_alpha = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value <= 3.41
            and self.u_value > 1.7
            and self.solar_heat_gain_coef <= 0.3
            and self.solar_heat_gain_coef > 0.25
        ):  # Zone 14
            Ts1 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)], axis=0)
            Rs1 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)], axis=0)
            Ts2 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs2 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.solar_heat_gain_coef - 0.25) / (
                0.3 - 0.25
            )
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.solar_heat_gain_coef - 0.25) / (
                0.3 - 0.25
            )
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value <= 3.41
            and self.u_value > 1.7
            and self.solar_heat_gain_coef <= 0.25
        ):  # Zone 15
            Ts_alpha = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)], axis=0)
            Rs_alpha = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)], axis=0)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef > 0.65
        ):  # Zone 16 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = Ts_A(Cosalpha)
            Rs2 = Rs_A(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef <= 0.65
            and self.solar_heat_gain_coef > 0.6
        ):  # Zone 17 bilinear interpolation
            Ts11 = Ts_E(Cosalpha)
            Rs11 = Rs_E(Cosalpha)
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs21 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts22 = Ts_A(Cosalpha)
            Rs22 = Rs_A(Cosalpha)
            Ts_alpha = (
                Ts11 * (0.65 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Ts21 * (self.solar_heat_gain_coef - 0.6) * (4.54 - self.u_value)
                + Ts12 * (0.65 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Ts22 * (self.solar_heat_gain_coef - 0.6) * (self.u_value - 3.41)
            ) / ((0.65 - 0.6) * (4.54 - 3.41))
            Rs_alpha = (
                Rs11 * (0.65 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Rs21 * (self.solar_heat_gain_coef - 0.6) * (4.54 - self.u_value)
                + Rs12 * (0.65 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Rs22 * (self.solar_heat_gain_coef - 0.6) * (self.u_value - 3.41)
            ) / ((0.65 - 0.6) * (4.54 - 3.41))
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef <= 0.6
            and self.solar_heat_gain_coef > 0.55
        ):  # Zone 18 linear interpolation
            Ts1 = Ts_E(Cosalpha)
            Rs1 = Rs_E(Cosalpha)
            Ts2 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs2 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef <= 0.55
            and self.solar_heat_gain_coef > 0.5
        ):  # Zone 19 bilinear interpolation
            Ts11 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs11 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts12 = Ts_E(Cosalpha)
            Rs12 = Rs_E(Cosalpha)
            Ts21 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs21 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts22 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs22 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts_alpha = (
                Ts11 * (0.55 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Ts21 * (self.solar_heat_gain_coef - 0.5) * (4.54 - self.u_value)
                + Ts12 * (0.55 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Ts22 * (self.solar_heat_gain_coef - 0.5) * (self.u_value - 3.41)
            ) / ((0.55 - 0.5) * (4.54 - 3.41))
            Rs_alpha = (
                Rs11 * (0.55 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Rs21 * (self.solar_heat_gain_coef - 0.5) * (4.54 - self.u_value)
                + Rs12 * (0.55 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Rs22 * (self.solar_heat_gain_coef - 0.5) * (self.u_value - 3.41)
            ) / ((0.55 - 0.5) * (4.54 - 3.41))
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef <= 0.5
            and self.solar_heat_gain_coef > 0.45
        ):  # Zone 20 linear interpolation
            Ts1 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs1 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts2 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs2 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef <= 0.45
            and self.solar_heat_gain_coef > 0.3
        ):  # Zone 21 bilinear interpolation
            Ts11 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs11 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts12 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs12 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts21 = Ts_D(Cosalpha)
            Rs21 = Rs_D(Cosalpha)
            Ts22 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs22 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts_alpha = (
                Ts11 * (0.45 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Ts21 * (self.solar_heat_gain_coef - 0.3) * (4.54 - self.u_value)
                + Ts12 * (0.45 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Ts22 * (self.solar_heat_gain_coef - 0.3) * (self.u_value - 3.41)
            ) / ((0.45 - 0.3) * (4.54 - 3.41))
            Rs_alpha = (
                Rs11 * (0.45 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Rs21 * (self.solar_heat_gain_coef - 0.3) * (4.54 - self.u_value)
                + Rs12 * (0.45 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Rs22 * (self.solar_heat_gain_coef - 0.3) * (self.u_value - 3.41)
            ) / ((0.45 - 0.3) * (4.54 - 3.41))
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef <= 0.3
            and self.solar_heat_gain_coef > 0.25
        ):  # Zone 22 bilinear interpolationm
            Ts11 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)], axis=0)
            Rs11 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)], axis=0)
            Ts12 = np.mean(
                [Ts_F(Cosalpha), Ts_G(Cosalpha), Ts_H(Cosalpha), Ts_I(Cosalpha)], axis=0
            )
            Rs12 = np.mean(
                [Rs_F(Cosalpha), Rs_G(Cosalpha), Rs_H(Cosalpha), Rs_I(Cosalpha)], axis=0
            )
            Ts21 = Ts_D(Cosalpha)
            Rs21 = Rs_D(Cosalpha)
            Ts22 = Ts_D(Cosalpha)
            Rs22 = Rs_D(Cosalpha)
            Ts_alpha = (
                Ts11 * (0.3 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Ts21 * (self.solar_heat_gain_coef - 0.25) * (4.54 - self.u_value)
                + Ts12 * (0.3 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Ts22 * (self.solar_heat_gain_coef - 0.25) * (self.u_value - 3.41)
            ) / ((0.3 - 0.25) * (4.54 - 3.41))
            Rs_alpha = (
                Rs11 * (0.3 - self.solar_heat_gain_coef) * (4.54 - self.u_value)
                + Rs21 * (self.solar_heat_gain_coef - 0.25) * (4.54 - self.u_value)
                + Rs12 * (0.3 - self.solar_heat_gain_coef) * (self.u_value - 3.41)
                + Rs22 * (self.solar_heat_gain_coef - 0.25) * (self.u_value - 3.41)
            ) / ((0.3 - 0.25) * (4.54 - 3.41))
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 3.41
            and self.u_value <= 4.54
            and self.solar_heat_gain_coef <= 0.25
        ):  # Zone 23 linear interpolation
            Ts1 = np.mean([Ts_F(Cosalpha), Ts_H(Cosalpha)], axis=0)
            Rs1 = np.mean([Rs_F(Cosalpha), Rs_H(Cosalpha)], axis=0)
            Ts2 = Ts_D(Cosalpha)
            Rs2 = Rs_D(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.u_value - 3.41) / (4.54 - 3.41)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif self.u_value > 4.5 and self.solar_heat_gain_coef > 0.65:  # Zone 24
            Ts_alpha = Ts_A(Cosalpha)
            Rs_alpha = Rs_A(Cosalpha)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 4.5
            and self.solar_heat_gain_coef > 0.6
            and self.solar_heat_gain_coef <= 0.65
        ):  # Zone 25 linear interpolation
            Ts1 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs1 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts2 = Ts_A(Cosalpha)
            Rs2 = Rs_A(Cosalpha)
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.solar_heat_gain_coef - 0.6) / (
                0.65 - 0.6
            )
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.solar_heat_gain_coef - 0.6) / (
                0.65 - 0.6
            )
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 4.5
            and self.solar_heat_gain_coef > 0.45
            and self.solar_heat_gain_coef <= 0.6
        ):  # Zone 26
            Ts_alpha = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs_alpha = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif (
            self.u_value > 4.5
            and self.solar_heat_gain_coef > 0.3
            and self.solar_heat_gain_coef <= 0.45
        ):  # Zone 27 linear interpolation
            Ts1 = Ts_D(Cosalpha)
            Rs1 = Rs_D(Cosalpha)
            Ts2 = np.mean(
                [Ts_B(Cosalpha), Ts_D(Cosalpha), Ts_C(Cosalpha), Ts_D(Cosalpha)], axis=0
            )
            Rs2 = np.mean(
                [Rs_B(Cosalpha), Rs_D(Cosalpha), Rs_C(Cosalpha), Rs_D(Cosalpha)], axis=0
            )
            Ts_alpha = Ts1 + (Ts2 - Ts1) * (self.solar_heat_gain_coef - 0.3) / (
                0.45 - 0.3
            )
            Rs_alpha = Rs1 + (Rs2 - Rs1) * (self.solar_heat_gain_coef - 0.3) / (
                0.45 - 0.3
            )
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        elif self.u_value > 4.5 and self.solar_heat_gain_coef <= 0.3:  # Zone 28
            Ts_alpha = Ts_D(Cosalpha)
            Rs_alpha = Rs_D(Cosalpha)
            Ts_abs_alpha = self.Ts * np.concatenate(([1], Ts_alpha))
            Rs_abs_alpha = self.Rs_f * np.concatenate(([1], Rs_alpha))
        else:
            print("Error in the Function. Check the Curve definition and areas")
        # Controlling the two vectors

        Ts_abs_alpha[9] = 0
        Rs_abs_alpha_2 = np.zeros(10)
        Rs_abs_alpha_2[0] = Rs_abs_alpha[0]
        if max(Rs_abs_alpha) > 1:
            for i in range(1, 10):
                Rs_abs_alpha_2[i] = Rs_abs_alpha[0] + (
                    Rs_abs_alpha[i] - Rs_abs_alpha[0]
                ) * (1 - Rs_abs_alpha[0]) / (max(Rs_abs_alpha) - Rs_abs_alpha[0])
            Rs_abs_alpha = Rs_abs_alpha_2
        else:
            Rs_abs_alpha[9] = 1
        for i in range(10):
            if Rs_abs_alpha[i] + Ts_abs_alpha[i] > 1:
                Rs_abs_alpha[i] = 1 - Ts_abs_alpha[i]
        self.Rs_abs_alpha = Rs_abs_alpha
        self.Ts_abs_alpha = Ts_abs_alpha
        self.As_abs_alpha = 1 - self.Rs_abs_alpha - self.Ts_abs_alpha
        self.solar_heat_gain_coef_abs_alpha = (
            self.Ts_abs_alpha + self.N * self.As_abs_alpha
        )
        self.alpha2 = np.linspace(0, 90, 10)
        self.solar_heat_gain_coef_profile = splrep(
            self.alpha2, self.solar_heat_gain_coef_abs_alpha, s=0
        )
        self.solar_heat_gain_coef_profile_conv = splrep(
            self.alpha2, self.N * self.As_abs_alpha, s=0
        )
        self.solar_heat_gain_coef_profile_rad = splrep(
            self.alpha2, self.Ts_abs_alpha, s=0
        )

    def __str__(self):
        return f"""
SimpleWindow: {self.name}
    U-value: {self.u_value} {units["U_value"]}
    Solar Heat Gain Coef.: {self.solar_heat_gain_coef} {units["solar_heat_gain_coefficient"]}"""
