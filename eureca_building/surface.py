"""
This module includes functions to model a 3D surface
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import numpy as np
import pyclipper as pc

from eureca_building.config import CONFIG
from eureca_building.construction import Construction
from eureca_building.window import SimpleWindow
from eureca_building.exceptions import (
    Non3ComponentsVertex,
    SurfaceWrongNumberOfVertices,
    WindowToWallRatioOutsideBoundaries,
    InvalidSurfaceType,
    NonPlanarSurface,
    NegativeSurfaceArea,
)
from eureca_building._geometry_auxiliary_functions import (
    check_complanarity,
    polygon_area,
    normal_versor_2,
    centroid,
    _project,
    _project_inv,
)


# %% Surface class
def delete_duplicates(lst):
    seen = set()
    for item in lst:
        # Convert the inner list to a tuple before checking for uniqueness
        item_tuple = tuple(item)
        seen.add(item_tuple)
    return tuple(seen)



class Surface:
    """Class surface checks the complanarity and calculates the area.
    Then calculates the azimuth and tilt of the surface and set a surface
    type depending on the tilt angle
    
    co planarity:
        https://www.geeksforgeeks.org/program-to-check-whether-4-points-in-a-3-d-plane-are-coplanar/

    the area is calculated from:
        https://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates    
    """

    __warning_azimuth_subdivisions = False
    __warning_height_subdivisions = False
    
    _discharge_coefficient_nat_vent = 0.6

    def __init__(
            self,
            name: str,
            vertices: tuple = ((0, 0, 0), (0, 0, 0), (0, 0, 0)),
            wwr=None,
            subdivisions_solar_calc=None,
            surface_type=None,
            construction=None,
            window=None,
            n_window_layers: int = 1
    ):
        """Creates the surface object. Checks all the inputs using properties setter methods

        Parameters
        ----------
        name : str
            Name.
        vertices : tuple, default ((0, 0, 0), (0, 0, 0), (0, 0, 0))
            List of vertices coordinates [m]. The default is ([0, 0, 0], [0, 0, 0], [0, 0, 0]).
        wwr : float, default None
            window to wall ratio (between  and 0 and 1). The default is 0.0.
        subdivisions_solar_calc : dict, default None
            Something like {
            'azimuth_subdivisions': 8,
            'height_subdivisions': 3,
            }
            keys:
            azimuth_subdivisions : int, optional
            Number of azimuth discretization for radiation purposes. The default is 8.
            height_subdivisions : int, optional
            Number of height discretization for radiation purposes. The default is 3.
        surface_type : str, default None
            Type of surface 'ExtWall' or 'GroundFloor' or 'Roof'.
            If not provided autocalculate.
        construction : eureca_building.construction.Construction
            the construction object with the materials
        window : eureca_building.window.SimpleWindow
            the Window object with the materials

        """

        self.name = name
        self._vertices = vertices
        self._centroid = centroid(self._vertices)

        # Area calculation

        self._area = polygon_area(self._vertices)

        # Considering only three points in calculating the normal vector could create
        # reverse orientations if the three points are in a non-convex angle of the surface
        #
        # for this reason theres an alternative way to calculate the normal,
        # implemented in function: normalAlternative
        #
        # reference: https://stackoverflow.com/questions/32274127/how-to-efficiently-determine-the-normal-to-a-polygon-in-3d-space


        self._normal = normal_versor_2(self._vertices)

        self._set_azimuth_and_zenith()

        if wwr is not None:
            self._wwr = wwr
        else:
            self._wwr = 0.0
        # Param Solar Calc
        if subdivisions_solar_calc is not None:
            self.subdivisions_solar_calc = subdivisions_solar_calc
        else:
            self.subdivisions_solar_calc = {"height_subdivisions": CONFIG.height_subdivisions,
                                            "azimuth_subdivisions": CONFIG.azimuth_subdivisions, }
        # Surfcae type
        if surface_type is None:
            self._set_auto_surface_type()
        else:
            self.surface_type = surface_type

        if construction is not None:
            self.construction = construction
        if window is not None:
            self.window = window
            
        # Window layout for natural ventilation
        self._define_windows_layout(n_window_layers=n_window_layers)

    @property
    def _vertices(self) -> tuple:
        return self.__vertices

    @_vertices.setter
    def _vertices(self, value: tuple):
        try:
            value = tuple(value)
        except ValueError:
            raise TypeError(f"Vertices of surface {self.name} are not a tuple: {value}")
        value = delete_duplicates(value)
        if len(value) < 3:  # Not a plane - no area
            raise SurfaceWrongNumberOfVertices(
                f"Surface {self.name}. Number of vertices lower than 3: {value}"
            )
        for vtx in value:
            if not isinstance(vtx, tuple):
                raise TypeError(
                    f"Vertices of surface {self.name} are not a tuple: {value}"
                )
            if len(vtx) != 3:
                raise Non3ComponentsVertex(
                    f"Surface {self.name} has a vertex with len() != 3: {value}"
                )
            try:
                float(vtx[0])
                float(vtx[1])
                float(vtx[2])
            except ValueError:
                raise ValueError(
                    f"Surface {self.name}. One vertex contains non float values: {vtx}"
                )
            # Check coplanarity

            if not check_complanarity(value):
                raise NonPlanarSurface(f"Surface {self.name}. Non planar points")
        self.__vertices = value

    @property
    def _area(self) -> float:
        return self.__area

    @_area.setter
    def _area(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Surface {self.name}, area is not an float: {value}")
        if value < 0.0:
            raise NegativeSurfaceArea(
                f"Surface {self.name}, negative surface area: {value}"
            )
        if float(value) == 0.0:
            self.__area = 1e-10
        else:
            self.__area = value

    @property
    def _wwr(self) -> float:
        return self.__wwr

    @_wwr.setter
    def _wwr(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Surface {self.name}, wwr is not an float: {value}")
        if value < 0.0 or value > 0.999:
            raise WindowToWallRatioOutsideBoundaries(
                f"Surface {self.name}, wwrS must included between 0 and 1: {value}"
            )
        self._calc_glazed_and_opaque_areas(value)
        self.__wwr = value

    @property
    def subdivisions_solar_calc(self) -> dict:
        return self._subdivisions_solar_calc

    @subdivisions_solar_calc.setter
    def subdivisions_solar_calc(self, value: dict):
        if not isinstance(value, dict):
            raise TypeError(
                f"Surface {self.name}, subdivisions_solar_calc must be a dict: {value}"
            )
        try:
            self._azimuth_subdivisions = value["azimuth_subdivisions"]
        except KeyError:
            raise KeyError(
                f"Surface {self.name}, subdivisions_solar_calc must contain an azimuth_subdivisions key: {value}"
            )
        try:
            self._height_subdivisions = value["height_subdivisions"]
        except KeyError:
            raise KeyError(
                f"Surface {self.name}, subdivisions_solar_calc must contain an height_subdivisions key: {value}"
            )
        self._subdivisions_solar_calc = value
        self._set_azimuth_and_zenith_solar_radiation()

    @property
    def _azimuth_subdivisions(self) -> int:
        return self.__azimuth_subdivisions

    @_azimuth_subdivisions.setter
    def _azimuth_subdivisions(self, value: int):
        try:
            value = int(value)
        except ValueError:
            raise TypeError(
                f"Surface {self.name}, azimuth_subdivisions is not an int: {value}"
            )
        if value < 1 or value > 100:
            # Check if unreasonable values provided
            raise ValueError(
                f"Surface {self.name}, azimuth_subdivisions must be > 1 and lower than 100: {value}"
            )
        if value > 16 and not self.__warning_azimuth_subdivisions:
            logging.warning(
                f"For one or more surfaces azimuth_subdivisions is high: {value}.\nThe calculation time can be long"
            )
            self.__warning_azimuth_subdivisions = True
        self.__azimuth_subdivisions = value

    @property
    def _height_subdivisions(self) -> int:
        return self.__height_subdivisions

    @_height_subdivisions.setter
    def _height_subdivisions(self, value: int):
        try:
            value = int(value)
        except ValueError:
            raise TypeError(
                f"Surface {self.name}, height_subdivisions is not an int: {value}"
            )
        if value < 1 or value > 50:
            # Check if unreasonable values provided
            raise ValueError(
                f"Surface {self.name}, height_subdivisions must be > 1 and lower than 50: {value}"
            )
        if value > 6 and not self.__warning_height_subdivisions:
            logging.warning(
                f"For one or more surfaces height_subdivisions is high: {value}.\nThe calculation time can be long"
            )
            self.__warning_height_subdivisions = True
        self.__height_subdivisions = value

    @property
    def surface_type(self):
        return self._surface_type

    @surface_type.setter
    def surface_type(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError(f"Surface {self.name}, surface_type is not a str: {value}")
        if value not in ["ExtWall", "GroundFloor", "Roof"]:
            raise InvalidSurfaceType(
                f"Surface {self.name}, surface_type must choosen from: [ExtWall, GroundFloor, Roof] {value}"
            )
        self._surface_type = value

    @property
    def construction(self):
        return self._construction

    @construction.setter
    def construction(self, value):
        if not isinstance(value, Construction):
            raise TypeError(f"Surface {self.name}, construction must be a Construction object: {type(value)}")
        self._construction = value

    @property
    def window(self):
        return self._window

    @window.setter
    def window(self, value):
        if not isinstance(value, SimpleWindow):
            raise TypeError(f"Surface {self.name}, window must be a SimpleWindow object: {type(value)}")
        self._window = value

    def _set_azimuth_and_zenith(self):
        """Internal method to calculate azimuth and zenith
        """

        # set the azimuth and zenith

        if self._normal[2] == 1:
            self._height = 0
            self._azimuth = 0
        elif self._normal[2] == -1:
            self._height = 180
            self._azimuth = 0
        else:
            self._height = 90 - np.degrees(
                np.arctan(
                    (
                            self._normal[2]
                            / (np.sqrt(self._normal[0] ** 2 + self._normal[1] ** 2))
                    )
                )
            )
            if self._normal[1] == 0:
                if self._normal[0] > 0:
                    self._azimuth = -90
                elif self._normal[0] < 0:
                    self._azimuth = 90
            else:
                if self._normal[1] < 0:
                    self._azimuth = np.degrees(
                        np.arctan(self._normal[0] / self._normal[1])
                    )
                else:
                    if self._normal[0] < 0:
                        self._azimuth = 180 + np.degrees(
                            np.arctan(self._normal[0] / self._normal[1])
                        )
                    else:
                        self._azimuth = -180 + np.degrees(
                            np.arctan(self._normal[0] / self._normal[1])
                        )

    def _calc_glazed_and_opaque_areas(self, wwr):
        """Internal method to calculate glazed and opaque ares

        Parameters
        ----------
        wwr : float
            Window-to-wall ration. Number between 0 and 1

        """
        self._opaque_area = (1 - wwr) * self._area
        self._glazed_area = wwr * self._area
        
    def _define_windows_layout(self, n_window_layers: int = 1):
        """USED FOR NATURAL VENTILATION
        Defines the windows layout (sill height, width, number, ...)

        Parameters
        ----------
        n_window_layers : int, default 1
            Number of rows to consider
        """
        _h_window_default = 1.5
        if not isinstance(n_window_layers, int):
            raise TypeError(f"Surface {self.name}: number of window layers is not an integer: n_window_layers {n_window_layers}")
        area_layer = self._glazed_area/n_window_layers
        self._w_window = area_layer/_h_window_default
        self._h_window = _h_window_default
        self._h_bottom_windows = np.array([1.2 + n*3.3 for n in range(n_window_layers)])       
        self._h_top_windows = self._h_bottom_windows + self._h_window
        self._a_coeff = self._discharge_coefficient_nat_vent*self._w_window   # Coeff a = discharge coeff * width window for calculation of natural ventilation flow rate

    def _set_azimuth_and_zenith_solar_radiation(self):
        """Internal method to calculate rounded azimuth and zenith
        """
        # Azimuth and tilt approximation

        delta_a = 360 / (2 * self._azimuth_subdivisions)
        delta_h = 90 / (2 * self._height_subdivisions)
        x = np.arange(-delta_h, 90 + 2 * delta_h, 2 * delta_h)

        for n in range(len(x) - 1):
            if self._height >= x[n] and self._height < x[n + 1]:
                self._height_round = int((x[n] + x[n + 1]) / 2)
                self._sky_view_factor = (1 + np.cos(np.radians(self._height_round))) / 2
            elif self._height >= x[-1] and self._height < 150:
                self._height_round = 90
                self._sky_view_factor = (1 + np.cos(np.radians(self._height_round))) / 2
            else:
                self._height_round = 0  # Only to avoid errors
                self._sky_view_factor=1 # Also to avoid errors                
        y = np.arange(-180 - delta_a, 180 + 2 * delta_a, 2 * delta_a)
        for n in range(len(y) - 1):
            if self._azimuth >= y[n] and self._azimuth < y[n + 1]:
                self._azimuth_round = int((y[n] + y[n + 1]) / 2)
                if self._azimuth_round == 180:
                    self._azimuth_round = -180
        if self._height_round == 0:
            self._azimuth_round = 0

    def _set_auto_surface_type(self):
        """Uses tilt to autoset surface type.
        tilt > 150 --> GroundFloor
        tilt < 40 --> Roof
        else --> ExtWall
        """
        # Set surface inclination

        if self._height < 40:
            self.surface_type = "Roof"
        elif self._height > 150:
            self.surface_type = "GroundFloor"
        else:
            self.surface_type = "ExtWall"

    def max_height(self):
        """Calculates max height from the most high vertex
        """
        hmax = 0
        for vert in self.__vertices:
            hmax = max(hmax, vert[2])
        return hmax

    def min_height(self):
        """Calculates max height from the most low vertex
        """
        hmin = 10000
        for vert in self.__vertices:
            hmin = min(hmin, vert[2])
        return hmin

    def get_VDI6007_surface_params(self, asim=None):
        """Calculates R and C using VDI6007 method.

        Parameters
        ----------
        asim : bool
            Whether the surface is asimmetric (True) or not (False)

        Returns
        -------
        tuple
            R, C -> Thermal Resistance and Capacity
        """
        if asim is None:
            if self.surface_type in ["ExtWall", "GroundFloor", "Roof"]:
                asim = True
            else:
                asim = False
        try:
            R1, C1 = self.construction._VDI6007_surface_params(self._area, asim)
        except AttributeError:
            raise AttributeError(
                f"Surface {self.name}, construction not specified"
            )
        return R1, C1

    def get_surface_external_radiative_coefficient(self):
        """Returns the radiative heat exchange coefficient.

        Returns
        -------
        float
        """
        # From standard average value
        return 5 * 0.9  # W/(m2 K)

    def check_surface_coincidence(self, other_surface):
        """Check if two surface are coincident returning True or False

        Parameters
        ----------
        other_surface : eureca_building.surface.Surface
            another surface object

        Returns
        -------
        bool
            Are the surfaces coincident? True/False
        """

        # Check Input data type

        if not isinstance(other_surface, Surface):
            raise ValueError(
                f"ERROR Surface class, surface {self.name}, check_surface_coincidence. other_surface is not a Surface object: otherSurface {other_surface}")

        # Check the coincidence of two surface looking firstly at the coplanarity
        # of the points and then the direction of the normals vectors

        flagPoints = False
        plane = list(self._vertices)

        # Coplanarity test
        for i in other_surface._vertices:
            if check_complanarity(plane + [i], precision=5):
                flagPoints = True

        # Normal vector test
        flagNormal = False
        if np.linalg.norm(self._normal + other_surface._normal) < 0.2:
            flagNormal = True
        return (flagNormal and flagPoints)

    def calculate_intersection_area(self, other_surface):
        '''Calculates the area between two adjacent surfaces
        reference: https://stackoverflow.com/questions/39003450/transform-3d-polygon-to-2d-perform-clipping-and-transform-back-to-3d

        Parameters
        ----------
        other_surface : eureca_building.surface.Surface
            another surface object

        Returns
        -------
        float
            The intersection area [m2]
        '''

        # Check Input data type

        if not isinstance(other_surface, Surface):
            raise ValueError(
                f"ERROR Surface class, surface {self.name}, calculate_intersection_area. other_surface is not a Surface object: otherSurface {other_surface}")

        # Check the coincidence of two surface looking firstly at the coplanarity
        # of the points and then the direction of the normals vectors
        a = self._normal[0] * self._vertices[0][0] + self._normal[1] * self._vertices[0][1] + self._normal[2] * \
            self._vertices[0][2]
        proj_axis = max(range(3), key=lambda i: abs(self._normal[i]))
        projA = [_project(x, proj_axis) for x in self._vertices]
        projB = [_project(x, proj_axis) for x in other_surface._vertices]
        scaledA = pc.scale_to_clipper(projA)
        scaledB = pc.scale_to_clipper(projB)
        clipper = pc.Pyclipper()
        clipper.AddPath(scaledA, poly_type=pc.PT_SUBJECT, closed=True)
        clipper.AddPath(scaledB, poly_type=pc.PT_CLIP, closed=True)
        intersections = clipper.Execute(pc.CT_INTERSECTION, pc.PFT_NONZERO, pc.PFT_NONZERO)
        intersections = [pc.scale_from_clipper(i) for i in intersections]
        if len(intersections) == 0:
            return 0
        intersection = tuple([_project_inv(x, proj_axis, a, self._normal) for x in intersections[0]])
        area = polygon_area(intersection)
        return area if area > 0 else 0.

    def reduce_area(self, area_to_reduce):
        '''Reduces the area of the surface by an input area

        Parameters
        ----------
        area_to_reduce : float
            the area to subtract [m2] to the total area
        '''

        # Check Input data type

        if not isinstance(area_to_reduce, float) or area_to_reduce < 0.:
            try:
                area_to_reduce = float(area_to_reduce)
            except ValueError:
                raise ValueError(
                    f"ERROR Surface class, surface {self.name}, reduce_area. The area is not a positive float: AreaToReduce {area_to_reduce}")

        # Area reduction

        if self._area - area_to_reduce > 0.0000001:
            self._area = self._area - area_to_reduce
        else:
            self._area = 0.0000001
        self._calc_glazed_and_opaque_areas(self._wwr) # This runs again the glazed and opaque calculation

    def __str__(self):
        return f"""
Name: {self.name}
    Type: {self.surface_type}
    Azimuth: {self._azimuth:.2f}
    Height: {self._height:.2f}
    U value: {self.construction._u_value:.2f}
    Area: {self._area} ({self._wwr:.1%} glazed)
"""

# %%---------------------------------------------------------------------------------------------------
# %% SurfaceInternalMass class


class SurfaceInternalMass:
    """Class to define a surface for thermal capacity using area and surface type
    with a specific geometry
    """

    def __init__(self, name: str, area: float = 0., surface_type=None, construction=None):
        """It creates the SurfaceInternalMass object, like the Surface class, but without vertexes and geometry

        Parameters
        ----------
        name : string
            name of the surface
        area : float, default 0.
            area of the internal surface
        surface_type : str, default None
            Type of internal surface: 'IntWall' or 'IntCeiling'
        construction : eureca_building.construction.Construction
            The construction to be assigned to the SurfaceInternalMass

        """

        # Check input data type
        self.name = name
        self._area = area
        self.surface_type = surface_type
        if construction is not None:
            self.construction = construction

    @property
    def _area(self) -> float:
        return self.__area

    @_area.setter
    def _area(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"SurfaceInternalMass {self.name}, area is not an float: {value}")
        if value < 0.0:
            raise NegativeSurfaceArea(
                f"SurfaceInternalMass {self.name}, negative surface area: {value}"
            )
        if float(value) == 0.0:
            self.__area = 1e-10
        else:
            self.__area = value
        self._opaque_area = self._area
        self._glazed_area = 0.

    @property
    def surface_type(self):
        return self._surface_type

    @surface_type.setter
    def surface_type(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError(f"SurfaceInternalMass {self.name}, surface_type is not a str: {value}")
        if value == None:
            logging.warning(f"SurfaceInternalMass {self.name}, surface_type is None: {value}. IntWall will be assigned")
            value = "IntWall"
        if value not in ["IntWall", "IntCeiling", "IntFloor"]:
            raise InvalidSurfaceType(
                f"SurfaceInternalMass {self.name}, surface_type must choosen from: [IntWall, IntCeiling, IntFloor] {value}"
            )
        self._surface_type = value

    @property
    def construction(self):
        return self._construction

    @construction.setter
    def construction(self, value):
        if not isinstance(value, Construction):
            raise TypeError(
                f"SurfaceInternalMass {self.name}, construction must be a Construction object: {type(value)}")
        self._construction = value

    def get_VDI6007_surface_params(self, asim=False):
        """Calculates R and C using VDI6007 method.

        Parameters
        ----------
        asim : bool
            Whether the surface is asimmetric (True) or not (False)

        Returns
        -------
        tuple
            R, C -> Thermal Resistance and Capacity
        """
        try:
            R1, C1 = self.construction._VDI6007_surface_params(self._opaque_area, asim)
        except AttributeError:
            raise AttributeError(
                f"SurfaceInternalMass {self.name}, construction not specified"
            )
        return R1, C1
