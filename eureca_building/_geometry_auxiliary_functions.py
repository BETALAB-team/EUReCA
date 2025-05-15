"""
Internal module with auxiliary functions for polygon geometry calculations.

Includes:
- Normal vector estimation
- Polygon area and centroid
- Coplanarity checks
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import numpy as np


def normal_versor(a, b, c):
    """This function starting from three points defines the normal vector of the plane
    
    Parameters
    ----------
    a : list
        list with three floats with the coordinates of the first point
    b : list
        list with three floats with the coordinates of the second point
    c : list
        list with three floats with the coordinates of the third point
 
    Returns
    -------
    tuple
        the three components of the normal vector (anticlockwise)

    Raises
    ------
    TypeError
        if input not Lists
    ValueError
        if number not floats
    """

    # Check input data type

    for ipt in (a, b, c):
        if not isinstance(ipt, list) and not isinstance(ipt, tuple):
            raise TypeError(
                f"ERROR normal_versor function, an input is not a list: input {ipt}"
            )
        if len(ipt) != 3:
            raise TypeError(
                f"ERROR normal_versor function, a vertex is not a list of 3 components: input {ipt}"
            )
        try:
            ipt = list(ipt)
            ipt[0] = float(ipt[0])
            ipt[1] = float(ipt[1])
            ipt[2] = float(ipt[2])
        except ValueError:
            raise ValueError(
                f"ERROR normal_versor function, a coordinate is not a float: input {ipt}"
            )
    # unit normal vector of plane defined by points a, b, and c
    x = np.linalg.det([[1, a[1], a[2]], [1, b[1], b[2]], [1, c[1], c[2]]])
    y = np.linalg.det([[a[0], 1, a[2]], [b[0], 1, b[2]], [c[0], 1, c[2]]])
    z = np.linalg.det([[a[0], a[1], 1], [b[0], b[1], 1], [c[0], c[1], 1]])
    magnitude = (x ** 2 + y ** 2 + z ** 2) ** 0.5
    return (x / magnitude, y / magnitude, z / magnitude)


# %%


def normal_versor_2(vert_list):
    """Alternative: This function starting from three points defines the normal vector of the plane
    
    Parameters
    ----------
    vert_list : tuple
        list with n lists of three floats (n vertices)
 
    Returns
    -------
    numpy.array
        coordinates of the centroid (3 components)

    Raises
    ------
    TypeError
        if input not tuples, and not long 3 components
    ValueError
        if number not floats
    """

    # Check input data type

    if not isinstance(vert_list, tuple):
        raise TypeError(
            f"ERROR normal_versor_2 function, the input is not a tuple: input {vert_list}"
        )
    # if len(vert_list) > 3:
    # wrn(f"normalAlternative function, there vertlist should be 3 components long: vertList {vertList}")

    for vtx in vert_list:
        if not isinstance(vtx, tuple):
            raise TypeError(
                f"ERROR normalAlternative function, an input is not a tuple: input {vtx}"
            )
        if len(vtx) != 3:
            raise TypeError(
                f"ERROR normalAlternative function, a vertex is not a tuple of 3 components: input {vtx}"
            )
        try:
            float(vtx[0])
            float(vtx[1])
            float(vtx[2])
        except ValueError:
            raise ValueError(
                f"ERROR normalAlternative function, a coordinate is not a float: input {vtx}"
            )
    c = centroid(vert_list)
    crossProd = np.array([0.0, 0, 0])
    for i in range(len(vert_list)):
        a = np.array(vert_list[i - 1]) - c
        b = np.array(vert_list[i]) - c
        crossProd += np.cross(a, b)
    return crossProd / np.linalg.norm(crossProd)





def polygon_area(poly):
    """From a list of points calculates the area of the polygon

    Parameters
    ----------
    poly: list
        list of list of floats (polygon) with three floats with the coordinates of the first point

    Returns
    -------
    float
        Area [m2]

    Raises
    ------
    TypeError
        if input not tuples, and not long 3 components
    ValueError
        if number not floats
    """

    # Check input data type

    if not isinstance(poly, tuple):
        raise TypeError(
            f"ERROR polygon_area function, the input is not a list: input {poly}"
        )
    for vtx in poly:
        if not isinstance(vtx, list) and not isinstance(vtx, tuple):
            raise TypeError(
                f"ERROR polygon_area function, an input is not a list: input {vtx}"
            )
        if len(vtx) != 3:
            raise TypeError(
                f"ERROR polygon_area function, a vertex is not a list of 3 components: input {vtx}"
            )
        try:
            vtx = list(vtx)
            vtx[0] = float(vtx[0])
            vtx[1] = float(vtx[1])
            vtx[2] = float(vtx[2])
        except ValueError:
            raise ValueError(
                f"ERROR unit_normal function, a coordinate is not a float: input {vtx}"
            )
    # area of polygon poly

    if len(poly) < 3:  # Not a plane - no area
        logging.error("WARNING number of vertices lower than 3, the area will be zero")
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i + 1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, normal_versor_2((poly[0], poly[1], poly[2])))
    return float(abs(result / 2))


# %%


def check_complanarity(vert_list_tot, precision=1):
    """Checks the co-planarity of a list of points
    
    Parameters
    ----------
    vert_list_tot : list
        list of list of floats (polygon). list with n lists of three floats (n vertices)
    precision : float
        defines the tolerance of the control [m]
 
    Returns
    -------
    bool
        are they in the same plane? True or False

    Raises
    ------
    ValueError
        if precision not float
    """

    # Check input data type

    # if not isinstance(vert_list_tot, list):
    #     raise TypeError(
    #         f"ERROR check_complanarity function, the input is not a list: input {vert_list_tot}"
    #     )
    # for vtx in vert_list_tot:
    #     if not isinstance(vtx, list):
    #         raise TypeError(
    #             f"ERROR check_complanarity function, an input is not a list: input {vtx}"
    #         )
    #     if len(vtx) != 3:
    #         raise TypeError(
    #             f"ERROR check_complanarity function, a vertex is not a list of 3 components: input {vtx}"
    #         )
    #     try:
    #         vtx[0] = float(vtx[0])
    #         vtx[1] = float(vtx[1])
    #         vtx[2] = float(vtx[2])
    #     except ValueError:
    #         raise ValueError(
    #             f"ERROR check_complanarity function, a coordinate is not a float: input {vtx}"
    #         )
    try:
        precision = float(precision)
    except ValueError:
        raise ValueError(
            f"ERROR check_complanarity function, precision is not a float: precision {precision}"
        )
    # Look if they are coplanar

    flag = True
    for i in range(len(vert_list_tot) - 3):
        vert_list = vert_list_tot[i: (i + 4)]
        a1 = vert_list[1][0] - vert_list[0][0]
        b1 = vert_list[1][1] - vert_list[0][1]
        c1 = vert_list[1][2] - vert_list[0][2]
        a2 = vert_list[2][0] - vert_list[0][0]
        b2 = vert_list[2][1] - vert_list[0][1]
        c2 = vert_list[2][2] - vert_list[0][2]
        a = b1 * c2 - b2 * c1
        b = a2 * c1 - a1 * c2
        c = a1 * b2 - b1 * a2
        d = -a * vert_list[0][0] - b * vert_list[0][1] - c * vert_list[0][2]

        # equation of plane is: a*x + b*y + c*z = 0
        # checking if the 4th point satisfies
        # the above equation
        if not (
                np.abs(a * vert_list[3][0] + b * vert_list[3][1] + c * vert_list[3][2] + d)
                < precision
        ):
            flag = False
    return flag


# %%


def centroid(vert_list):
    """From a list of points calculates the centroid
    
    Parameters
    ----------
    vert_list : list
        list with n lists of three floats (n vertices)
 
    Returns
    -------
    numpy.array
        coordinates of the centroid (3 components)
    
    """

    # Check input data type

    if not isinstance(vert_list, tuple):
        raise TypeError(
            f"ERROR centroid function, the input is not a list: input {vert_list}"
        )
    for vtx in vert_list:
        if not isinstance(vtx, tuple):
            raise TypeError(
                f"ERROR centroid function, an input is not a list: input {vtx}"
            )
        if len(vtx) != 3:
            raise TypeError(
                f"ERROR centroid function, a vertex is not a list of 3 components: input {vtx}"
            )
        try:
            float(vtx[0])
            float(vtx[1])
            float(vtx[2])
        except ValueError:
            raise ValueError(
                f"ERROR centroid function, a coordinate is not a float: input {vtx}"
            )
    # Centroid calculation
    c = np.array([0.0, 0, 0], dtype=float)
    for i in vert_list:
        c += np.array(i)
    return c / len(vert_list)


# %%


def _project(x, proj_axis):
    """Internal Function: Project onto either the xy, yz, or xz plane. (We choose the one that avoids degenerate configurations, which is the purpose of proj_axis.)
    # In this example, we would be projecting onto the xz plane.
    https://stackoverflow.com/questions/39003450/transform-3d-polygon-to-2d-perform-clipping-and-transform-back-to-3d


    Parameters
    ----------
    x : list
        list of the three coordinates of a point
    proj_axis : int
        0, 1, or 2, is the axis where to project (x, y, z)

    Returns
    -------
    tuple
        tuple of three float, projected point

    """
    return tuple(c for i, c in enumerate(x) if i != proj_axis)


# %%


def _project_inv(x, proj_axis, a, v):
    """Internal Function: Returns the vector w in the walls' plane such that project(w) equals x.
    https://stackoverflow.com/questions/39003450/transform-3d-polygon-to-2d-perform-clipping-and-transform-back-to-3d


    Parameters
    ----------
    x : list
        list of the three coordinates of a point
    proj_axis : int
        0, 1, or 2, is the axis where to project (x, y, z)
    a : float
        a number so that v[0] * x + v[1] * y + v[2] * z = a is the equation of the plane containing your walls
    v : list
        list of three floats. the normal versor to the surface

    Returns
    -------
    tuple
        tuple of three float, projected point

    """
    w = list(x)
    w[proj_axis:proj_axis] = [0.0]
    c = a
    for i in range(3):
        c -= w[i] * v[i]
    c /= v[proj_axis]
    w[proj_axis] = c
    return tuple(w)
