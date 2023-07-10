"""
List of custom exceptions
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


class PropertyOutsideBoundaries(Exception):
    """
    Class to raise the exception property outside boundaries
    """

    def __init__(self, mat, prop, lim=None, unit=None, value=None):
        """Create an exception for properties problems

        Parameters
        ----------
        mat : str
            name of the property
        prop : str
            string with the property
        unit : str
            unit of the property
        lim : list
            list with boundaries
        """

        super().__init__(mat, prop, lim, unit, value)
        self.prop = prop
        self.mat = mat
        self.lim = lim
        self.unit = unit
        self.value = value


class MaterialPropertyOutsideBoundaries(PropertyOutsideBoundaries):
    """
    Class to raise the exception material property outside boundaries
    """

    pass


class MaterialPropertyNotFound(Exception):
    """
    Class to raise the exception material property not found
    """

    pass


class WrongConstructionType(Exception):
    """
    Class to raise the exception wrong construction type
    """

    pass


class WrongMaterialType(Exception):
    """
    Class to raise the exception wrong construction type
    """

    pass


# %% Surfaces exeptions
class Non3ComponentsVertex(Exception):
    pass


class SurfaceWrongNumberOfVertices(Exception):
    pass


class WindowToWallRatioOutsideBoundaries(Exception):
    pass


class InvalidSurfaceType(Exception):
    pass


class NonPlanarSurface(Exception):
    pass


class NegativeSurfaceArea(Exception):
    pass


class InvalidScheduleType(Exception):
    pass


class ScheduleOutsideBoundaryCondition(Exception):
    pass


class InvalidScheduleDimension(Exception):
    pass


class ScheduleLengthNotConsistent(Exception):
    pass


class ConvectiveRadiantFractionError(Exception):
    pass


class InvalidHeatGainUnit(Exception):
    pass


class InvalidHeatGainSchedule(Exception):
    pass


class AreaNotProvided(Exception):
    pass


class PeopleNotProvided(Exception):
    pass

class SetpointTypeNotAllowed(Exception):
    pass

class SimulationError(Exception):
    pass
