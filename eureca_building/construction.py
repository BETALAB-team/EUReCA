"""
This module includes class for the typical construction
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import pandas as pd
import numpy as np

from eureca_building.exceptions import (
    MaterialPropertyOutsideBoundaries,
    WrongConstructionType,
)
from eureca_building.material import Material, AirGapMaterial
from eureca_building.units import units


# %% OpaqueMaterial class


class Construction(object):
    """
        Represents a multi-layer building construction, including methods to calculate:
            - Thermal resistance and U-values
            - Thermal capacitance (ISO 13790, VDI 6007)
            - Transfer matrix parameters for dynamic simulation
            
        Used for walls, roofs, floors, and internal partitions.
    """
    # Class attributes

    tot_heat_trans_coef = pd.DataFrame(
        {"Outside": [25, 25, 1000, 7.7, 6.7, 6.7], "Inside": [7.7, 7.7, 7.7, 7.7, 6.7, 6.7]},
        index=["ExtWall", "Roof", "GroundFloor", "IntWall", "IntCeiling", "IntFloor"],
    )

    rad_heat_trans_coef = 5.0

    def __init__(
            self, name: str, materials_list: list, construction_type: str = "ExtWall"
    ):
        """Initializes the Construction object

        Parameters
        ----------
        name : string
            name
        materials_list : list
            list of Materials or AirGapMaterials objects (Outside -> Inside)
        construction_type : string
            Choose from ["ExtWall", "Roof", "GroundFloor", "IntWall", "IntCeiling"]

        """

        # Check input data type
        for mat in materials_list:
            if (not isinstance(mat, Material)) and (
                    not isinstance(mat, AirGapMaterial)
            ):
                raise TypeError(
                    f"Construction {name}. materials_list must be a list of Materials or AirGapMaterial objects. Material {mat.name}"
                )
        if construction_type not in [
            "ExtWall",
            "Roof",
            "GroundFloor",
            "IntWall",
            "IntCeiling",
            "IntFloor",
        ]:
            raise WrongConstructionType(
                f'Construction {name}. construction type {construction_type} not in ["ExtWall", "Roof", "GroundFloor", "IntWall", "IntCeiling", "IntFloor"]'
            )
        self.name = name
        self.construction_type = construction_type
        self.materials_list = materials_list

        # Set outside and inside convective and radiant heat transfer coefficients

        self._R_si = 1 / (
            self.tot_heat_trans_coef.loc[self.construction_type]["Inside"]
        )
        self._R_se = 1 / (
            self.tot_heat_trans_coef.loc[self.construction_type]["Outside"]
        )
        self._conv_heat_trans_coef_int = (
                self.tot_heat_trans_coef.loc[self.construction_type]["Inside"]
                - self.rad_heat_trans_coef
        )
        self._conv_heat_trans_coef_ext = (
                self.tot_heat_trans_coef.loc[self.construction_type]["Outside"]
                - self.rad_heat_trans_coef
        )

        self.ext_absorptance = self.materials_list[0].thermal_absorptance
        self.number_of_layers = len(self.materials_list)

        # Calculation of the U-values [W/(m2 K)] and creation of some arrays
        self.net_thermal_resistance = 0
        self.thicknesses = np.zeros(self.number_of_layers)
        self.conductivities = np.zeros(self.number_of_layers)
        self.densities = np.zeros(self.number_of_layers)
        self.spec_heats = np.zeros(self.number_of_layers)
        self.thermal_resistances = np.zeros(self.number_of_layers)
        for i, mat in zip(range(self.number_of_layers), self.materials_list):
            self.net_thermal_resistance += mat.thermal_resistance
            self.thicknesses[i] = mat.thick
            self.conductivities[i] = mat.cond
            self.densities[i] = mat.dens
            self.spec_heats[i] = mat.spec_heat
            self.thermal_resistances[i] = mat.thermal_resistance
        self.thermal_resistance = self.net_thermal_resistance + self._R_si + self._R_se

        self._u_value_net = 1 / self.net_thermal_resistance
        self._u_value = 1 / self.thermal_resistance

        if self.construction_type == "GroundFloor":
            self._u_value = self._u_value * 0.7
        # Run ISO13790params and vdi6007params to calculate further parameters

        self._ISO13790_params()
        self._VDI6007_params()

    def _ISO13790_params(self):
        """Calculates ISO13790 params: k_int, k_est
        """

        # Set some parameters from the standard

        T = 86400
        sigma_2 = (T / np.pi) * (
                self.conductivities / (self.densities * self.spec_heats)
        )
        # Depth of penetration
        sigma = np.sqrt(sigma_2)
        eps = self.thicknesses / sigma

        # Thermal transfer matrix

        Z = np.zeros((2, 2, self.number_of_layers), complex)

        for i in range(self.number_of_layers):
            Z[0, 0, i] = np.cosh(eps[i]) * np.cos(eps[i]) + 1j * np.sinh(
                eps[i]
            ) * np.sin(eps[i])
            Z[1, 1, i] = Z[0, 0, i]
            Z[0, 1, i] = -(sigma[i] / (2 * self.conductivities[i])) * (
                    np.sinh(eps[i]) * np.cos(eps[i])
                    + np.cosh(eps[i]) * np.sin(eps[i])
                    + 1j
                    * (np.cosh(eps[i]) * np.sin(eps[i]) - np.sinh(eps[i]) * np.cos(eps[i]))
            )
            Z[1, 0, i] = -(self.conductivities[i] / sigma[i]) * (
                    np.sinh(eps[i]) * np.cos(eps[i])
                    - np.cosh(eps[i]) * np.sin(eps[i])
                    + 1j
                    * (np.sinh(eps[i]) * np.cos(eps[i]) + np.cosh(eps[i]) * np.sin(eps[i]))
            )
        Z_si = np.eye(2)
        # Internal surface resistance (convection and radiation, ISO 6946)
        Z_si[0, 1] = -self._R_si
        Z_se = np.eye(2)
        # External surface resistance (convection and radiation, ISO 6946)
        Z_se[0, 1] = -self._R_se
        i = self.number_of_layers - 1

        while i > 0:
            Z[:, :, i - 1] = np.matmul(Z[:, :, i], Z[:, :, i - 1])
            i = i - 1
        Z = Z[:, :, 0]
        Z = np.matmul(Z_se, Z)
        # Thermal transfer matrix
        Z = np.matmul(Z, Z_si)

        # Thermal capacitance
        self.k_est = (T / (2 * np.pi)) * np.abs((Z[0, 0] - 1) / Z[0, 1])
        # Thermal capacitance
        self.k_int = (T / (2 * np.pi)) * np.abs((Z[1, 1] - 1) / Z[0, 1])
        self.k_mean = (self.k_int + self.k_est) / 2

    def _VDI6007_params(self):
        """Calculates vdi6007 params

        Section 6.3

        vdi6007params Calculates the parameters (thermal resistance and thermal
        capacitance) associated with the building envelope LP according to the 2C
        model of standard VDI 6007

        Inputs (already attributes of the class)
          Matrix LP describes building envelope; each row is a building element.
          Columns are thickness (s), thermal conductivity (cond), density (rho)
          and specific heat (cp)
          Total surface area S of wall with building envelope LP
          Flag 'asim' indicates whether building component LP is asimmetrically
          loaded (asim = 1) or not (asim = 0) because in the first case C1_korr
          must be considered instead of C1

        Subscripts
        AW external walls and internal walls facing unheated areas
        IW internal walls

        Period T_bt
        T_bt = 7 days for a single building component;
        T_bt = 2 days for building components where thermal storage masses are
        thermally covered on the room side (eg suspended ceilings)
        Calculation are conducted with both periods, then the resulting
        parameters R1 and C1 are compared and the right T_bt is chosen according
        to the criterion given in the standard

        Outputs
           R1 - dynamic rhermal resistance [K/W]
           C1 - dynamic thermal capacity [K/W]
           Rw - static specific thermal resistance [(m2 K) / W]

        Determine thermal resistance R and thermal capacitance of layers, which are memorized as attributes
        """

        R = self.thermal_resistances  # layers thermal resistance [m2 K / W]
        # layers thermal capacitance [J m2 / K]
        C = self.spec_heats * self.densities * self.thicknesses

        T_bt = np.array([2, 7])  # days
        self.omega_bt = 2 * np.pi / (86400 * T_bt)
        # T_ra = 5        # days
        # omega_ra = 2*pi./(86400*T_ra)

        Z_t2 = np.zeros((2, 2, self.number_of_layers), complex)
        Z_t7 = np.zeros((2, 2, self.number_of_layers), complex)

        for i in range(self.number_of_layers):

            r = R[i]
            c = C[i]
            """
            R1 = np.zeros((2,1))
            R2 = R1
            C1 = R1
            C2 = R1 
            """
            # r_ah = self.massless[i]
            Av = np.zeros((2, 2, 2), complex)

            # if r_ah == 0:

            for om in range(2):
                arg = np.sqrt(0.5 * self.omega_bt[om] * r * c)

                Re_a11 = np.cosh(arg) * np.cos(arg)
                Im_a11 = np.sinh(arg) * np.sin(arg)
                Re_a22 = Re_a11
                Im_a22 = Im_a11

                Re_a12 = (
                        r
                        / (2 * arg)
                        * (np.cosh(arg) * np.sin(arg) + np.sinh(arg) * np.cos(arg))
                )
                Im_a12 = (
                        r
                        / (2 * arg)
                        * (np.cosh(arg) * np.sin(arg) - np.sinh(arg) * np.cos(arg))
                )

                Re_a21 = (
                        -arg / r * (np.cosh(arg) * np.sin(arg) - np.sinh(arg) * np.cos(arg))
                )
                Im_a21 = (
                        arg / r * (np.cosh(arg) * np.sin(arg) + np.sinh(arg) * np.cos(arg))
                )

                Av[0, 0, om] = Re_a11 + 1j * Im_a11
                Av[1, 1, om] = Re_a22 + 1j * Im_a22
                Av[0, 1, om] = Re_a12 + 1j * Im_a12
                Av[1, 0, om] = Re_a21 + 1j * Im_a21
            Z_t2[:, :, i] = Av[:, :, 0]
            Z_t7[:, :, i] = Av[:, :, 1]
        self._A1n_t2 = np.zeros((2, 2, 1), complex)
        self._A1n_t7 = np.zeros((2, 2, 1), complex)
        # From inside to outside
        # self._A1n_t2 = Z_t2[:, :, -1]
        # self._A1n_t7 = Z_t7[:, :, -1]
        # for t in range(-2, -self.number_of_layers - 1, -1):
        #     self._A1n_t2 = np.matmul(self._A1n_t2, Z_t2[:, :, t])
        #     self._A1n_t7 = np.matmul(self._A1n_t7, Z_t7[:, :, t])
        # From outside to inside
        self._A1n_t2 = Z_t2[:, :, 0]
        self._A1n_t7 = Z_t7[:, :, 0]
        for t in range(1, self.number_of_layers, 1):
            self._A1n_t2 = np.matmul(self._A1n_t2, Z_t2[:, :, t])
            self._A1n_t7 = np.matmul(self._A1n_t7, Z_t7[:, :, t])

    def _VDI6007_surface_params(self, sup, asim):
        """Calculates vdi6007 params, those which are calculated using the area of the surface

        Section 6.4

        Parameters
        ----------
        sup : float
            area of the surface [m2]
        asim : bool
            Is the surface non-adiabatic? True/False

        Returns
        -------
        tuple
            R1,C1: tuple of floats
            Resistance and capacitance R1 and C1
        """

        # Check input data type

        if not isinstance(sup, float):
            raise TypeError(
                f"ERROR surface {self.name} input sup is not a float: sup {sup}"
            )
        if not isinstance(asim, bool):
            raise TypeError(
                f"ERROR surface {self.name} input asim is not a boolean: asim {asim}"
            )
        # Procudeure Section 6.4

        if sup == 0:
            sup = 0.0000001
        rw = sum(self.thermal_resistances) / sup

        R1_t = dict()
        C1_t = dict()

        for a, omega, days in zip(
                [self._A1n_t2, self._A1n_t7],
                [self.omega_bt[0], self.omega_bt[1]],
                ["2", "7"],
        ):
            # rcValues Given the complex matrix of the building element BT, the function
            # calculates the values R1 and C1
            R1 = (
                    1
                    / sup
                    * (
                            (np.real(a[1, 1]) - 1) * np.real(a[0, 1])
                            + np.imag(a[1, 1]) * np.imag(a[0, 1])
                    )
                    / ((np.real(a[1, 1]) - 1) ** 2 + (np.imag(a[1, 1])) ** 2)
            )

            if asim == False:
                C1 = (
                        sup
                        * ((np.real(a[1, 1]) - 1) ** 2 + (np.imag(a[1, 1])) ** 2)
                        / (
                                omega
                                * (
                                        np.real(a[0, 1]) * np.imag(a[1, 1])
                                        - (np.real(a[1, 1]) - 1) * np.imag(a[0, 1])
                                )
                        )
                )
            else:
                # sarebbe C1_korr per pareti caricate asimmetricamente (pareti AW)
                C1 = (
                        sup
                        * (1 / (omega * R1 * sup))
                        * (
                                rw * sup
                                - np.real(a[0, 1]) * np.real(a[1, 1])
                                - np.imag(a[1, 1]) * np.imag(a[0, 1])
                        )
                        / (
                                np.real(a[1, 1]) * np.imag(a[0, 1])
                                - np.real(a[0, 1]) * np.imag(a[1, 1])
                        )
                )
            R1_t[days] = R1
            C1_t[days] = C1
        rr = R1_t["2"] / R1_t["7"]
        cr = C1_t["2"] / C1_t["7"]

        # versione di jacopo
        # if (rr>0.99 and cr<0.95) or (((rr<0.99 and cr<0.95) and (abs(rr-cr)>0.3))):
        if (rr > 0.99 and cr < 0.95) or (
                ((rr < 0.95 and cr < 0.95) and (abs(rr - cr) > 0.3))
        ):
            R1 = R1_t["2"]  # T_bt = 2 days
            C1 = C1_t["2"]
        else:
            R1 = R1_t["7"]  # T_bt = 7 days
            C1 = C1_t["7"]
        return R1, C1

    def __str__(self):
        """Just a print method
        """
        return f"""
Construction: {self.name}
    construction type: {self.construction_type}
    U-value: {self.U:.2f} {units["U_value"]}
    number of layers: {len(self.materials_list)}
    materias: {[str(mat.name) for mat in self.materials_list]}"""

    @classmethod
    def from_U_value(cls, name:str, u_value: float, weight_class:str = "Medium", construction_type:str ="ExtWall"):
        """This is a class method to create Construction object just from the U-value and weight class
        It creates just an equivalent material to reach the U-value

        For specific heat and density, the following assumptions are considered

        According to A.2.3 ISO 13786
        Mass class  Am [m²]	Cm [J/K]	k [J/(m² K)]	Depth penetration [m]	Spc heat [J/kg K]	rho [kg/m3]
        Very light	2.5	80000		32000		    0.1			            1000			    453
        Light		2.5	110000		44000		    0.1			            1000			    622
        Medium		2.5	165000		66000		    0.1			            1000			    933
        Heavy		3	260000		86666.66667	    0.1			            1000			    1226
        Very heavy	3.5	370000		105714.2857	    0.1			            1000	            1495

        Parameters
        ----------
        name : string
            name
        u_value : float
            u value of the construction [W/(m2 K)]
        weight_class : str
            class of weight from the following list: ["Very heavy", "Heavy, "Medium, "Light", "Very light"]
        construction_type : string
            Choose from ["ExtWall", "Roof", "GroundFloor", "IntWall", "IntCeiling"]

        Returns
        ----------
        eureca_building.construction.Construction
            Construction object from these values
            
        Notes
        -----
        This method creates an "equivalent material" with default absorptance and thickness (30 cm)
        based on ISO 13786 guidelines.

        """

        # Hypothesis 30 cm
        thickness = 0.3
        outside_ht_coef = cls.tot_heat_trans_coef.loc[construction_type]["Outside"]
        inside_ht_coef = cls.tot_heat_trans_coef.loc[construction_type]["Inside"]
        resistance = 1/u_value - 1/outside_ht_coef - 1/inside_ht_coef
        conductivity = 0.3/resistance

        """
	According to A.2.3 ISO 13786					
		Am [m²]	Cm [J/K]	k [J/(m² K)]	Depth penetration [m]	Spc heat [J/kg K]	rho [kg/m3]
Very light	2.5	80000		32000		    0.1			            1000			    453
Light		2.5	110000		44000		    0.1			            1000			    622
Medium		2.5	165000		66000		    0.1			            1000			    933
Heavy		3	260000		86666.66667	    0.1			            1000			    1226
Very heavy	3.5	370000		105714.2857	    0.1			            1000	            1495           
        """

        density = {
            "Very heavy": 1495.,
            "Heavy": 1226.,
            "Medium": 933.,
            "Light": 622.,
            "Very light": 453.,
        }[weight_class]
        equivalent_material = Material(
            f"Equivalent material construction {name}",
            thick=0.3,
            cond=conductivity,
            dens=density,
            thermal_absorptance=0.9,
        )

        return cls(name, [equivalent_material], construction_type)

