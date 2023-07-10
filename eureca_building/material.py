"""
This module includes classes and functions to implement materials properties
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

from eureca_building.exceptions import MaterialPropertyOutsideBoundaries
from eureca_building.units import units, material_limits


class Material:
    """Class used to define the material. Creates the material and checks the properties
    """

    name: str
    thick: float = 0.100  # Thickness [m]
    cond: float = 1.00  # Conductivity [W/mK]
    spec_heat: float = 1000.0  # Specific heat [J/kgK]
    dens: float = 1000.0  # Density [kg/m3]

    # Just to use the @property decorator and the setter function
    # _name: str = field(init = False, repr = False)
    # _thick: float = field(init = False, repr = False)
    # _cond: float = field(init = False, repr = False)
    # _spec_heat: float = field(init = False, repr = False)
    # _dens: float = field(init = False, repr = False)

    def __init__(
        self,
        name: str,
        thick: float = 0.100,
        cond: float = 1.00,
        spec_heat: float = 1000.0,
        dens: float = 1000.0,
        thermal_absorptance: float = 0.9,
    ):
        """Define the material and check the properties, using setter methods

        Parameters
        ----------
        name : str
            name
        thick : float, default 0.1
            thickness
        cond : float, default 1.
            conductivity
        spec_heat : float, default 1000.
            spec_heat
        dens : float, default 1000.
            density
        thermal_absorptance : float, default 0.9
            thermal absorptance  [-]

        Raises
        -------
        MaterialPropertyOutsideBoundaries
            If a material parameter is not allowed.
        """
        self.name = name
        self.thick = thick
        self.dens = dens
        self.cond = cond
        self.spec_heat = spec_heat
        self.thermal_absorptance = thermal_absorptance
        self.calc_params()

    @property
    def thick(self) -> float:
        return self._thick

    @thick.setter
    def thick(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Material {self.name}, thickness is not a float: {value}")
        if (
            value < material_limits["thickness"][0]
            or value > material_limits["thickness"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "thickness",
                lim=material_limits["thickness"],
                unit=units["length"],
                value=value,
            )
        self._thick = value

    @property
    def dens(self) -> float:
        return self._dens

    @dens.setter
    def dens(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Material {self.name}, density is not a float: {value}")
        if (
            value < material_limits["density"][0]
            or value > material_limits["density"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "density",
                lim=material_limits["density"],
                unit=units["density"],
                value=value,
            )
        self._dens = value

    @property
    def cond(self) -> float:
        return self._cond

    @cond.setter
    def cond(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, conductivity is not a float: {value}"
            )
        if (
            value < material_limits["conductivity"][0]
            or value > material_limits["conductivity"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "conductivity",
                lim=material_limits["conductivity"],
                unit=units["conductivity"],
                value=value,
            )
        self._cond = value

    @property
    def spec_heat(self) -> float:
        return self._spec_heat

    @spec_heat.setter
    def spec_heat(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, specific heat is not a float: {value}"
            )
        if (
            value < material_limits["specific_heat"][0]
            or value > material_limits["specific_heat"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "specific_heat",
                lim=material_limits["specific_heat"],
                unit=units["specific_heat"],
                value=value,
            )
        self._spec_heat = value

    @property
    def thermal_absorptance(self) -> float:
        return self._thermal_absorptance

    @thermal_absorptance.setter
    def thermal_absorptance(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, thermal_absorptance is not a float: {value}"
            )
        if (
            value < material_limits["absorptance"][0]
            or value > material_limits["absorptance"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "thermal_absorptance",
                lim=material_limits["absorptance"],
                unit=units["absorptance"],
                value=value,
            )
        self._thermal_absorptance = value

    def calc_capacity(self):
        """Thermal capacity calculation (thickness * density * specific_heat)
        """
        self.capacity = self.thick * self.dens * self.spec_heat

    def calc_resistance(self):
        """Thermal capacity calculation (thickness / conductivity)
        """
        self.thermal_resistance = self.thick / self.cond

    def calc_params(self):
        """Runs self.calc_capacity() and self.calc_resistance()
        """
        self.calc_capacity()
        self.calc_resistance()

    def __str__(self):
        """Return a beautiful parsed str of the properties

        Returns
        -------
        str
            A string to print properties
        """
        return f"""
Material: {self.name}
    thickness: {self.thick} {units["length"]}
    density: {self.dens} {units["density"]}
    specific heat: {self.spec_heat} {units["specific_heat"]}
    conductivity: {self.cond} {units["conductivity"]}"""


class AirGapMaterial:
    """A class used to define the air gap material
<    """

    name: str
    thick: float = 0.100  # Thickness [m]
    thermal_resistance: float = 1.00  # thermal_resistance [m2K/W]

    # Just to use the @property decorator and the setter function
    # _name: str = field(init = False, repr = False)
    # _thick: float = field(init = False, repr = False)
    # _resistance: float = field(init = False, repr = False)

    def __init__(
        self, name: str, thick: float = 0.100, thermal_resistance: float = 1.00
    ):
        """Defines the material and check the properties. Checks properties values using setter methods

        Parameters
        ----------
        name : str
            name
        thick : float, default 0.1
            thickness
        thermal_resistance : float, default 1.
            thermal_resistance

        Raises
        -------
        MaterialPropertyOutsideBoundaries
            If a material parameter is not allowed.
        """
        self.name = name
        self.thick = thick
        self.thermal_resistance = thermal_resistance

        # Just some equivalent values
        self.cond = self.thick / self.thermal_resistance
        self.dens = 1.2  # [kg/m3]
        self.spec_heat = 1005.0  # [J/kgK]

    @property
    def thick(self) -> float:
        return self._thick

    @thick.setter
    def thick(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Material {self.name}, thickness is not a float: {value}")
        if (
            value < material_limits["thickness"][0]
            or value > material_limits["thickness"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "thickness",
                lim=material_limits["thickness"],
                unit=units["length"],
                value=value,
            )
        self._thick = value

    @property
    def thermal_resistance(self) -> float:
        return self._thermal_resistance

    @thermal_resistance.setter
    def thermal_resistance(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Material {self.name}, thermal_resistance is not a float: {value}"
            )
        if (
            value < material_limits["thermal_resistance"][0]
            or value > material_limits["thermal_resistance"][1]
        ):
            # Value in [m]. Take a look to units
            # Check if thickenss is outside
            raise MaterialPropertyOutsideBoundaries(
                self.name,
                "thermal_resistance",
                lim=material_limits["thermal_resistance"],
                unit=units["thermal_resistance"],
                value=value,
            )
        self._thermal_resistance = value

    def __str__(self):
        """Return a beautiful parsed str of the properties

        Returns
        -------
        str
            A string to print properties
        """
        return f"""
AirGapMaterial: {self.name}
    thickness: {self.thick} {units["length"]}
    thermal resistance: {self.thermal_resistance} {units["thermal_resistance"]}"""
