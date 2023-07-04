"""
This module includes functions to manage the 2C model from the Standard VDI 6007
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import numpy as np


# %% ---------------------------------------------------------------------------------------------------
# %% Useful functions for VDI 6007 calculations

def impedence_parallel(R, C, T_RA=5.):
    '''Given two vectors (thermal resistances R and thermal
    capacitances C) of length m (number of walls of the same type, ie either
    IW or AW), calculates the equivalent complex thermal resistance Zeq
    according to T_RA (period in days)

    Parameters
    ----------
    R : numpy.array
        numpy.array with the surfaces resistances
    C: numpy.array
        numpy.array with the surfaces capacitances
    T_RA: float
        Reference time (number of days)

    Returns
    -------
    tuple
        tuple of floats: the equivalent resistance anc capacitance

    Raises
    ------
    TypeError
        if not numpy arrays or floats
    '''

    # Check input data type

    if not isinstance(R, np.ndarray):
        raise TypeError(f'ERROR impedenceParallel function, input R is not a np.array: R {R}')
    if not isinstance(C, np.ndarray):
        raise TypeError(f'ERROR impedenceParallel function, input C is not a np.array: C {C}')
    if not isinstance(T_RA, float):
        raise TypeError(f'ERROR impedenceParallel function, input T_RA is not a float: T_RA {T_RA}')

        # T_RA = 5;       % Period = 5 days
    omega_RA = 2 * np.pi / (86400 * T_RA)
    z = np.zeros(len(R), complex)
    z = (R + 1j / (omega_RA * C))  # vettore delle Z

    Z1eq = 1 / (sum(1 / z))  # Z equivalente

    R1eq = np.real(Z1eq)
    C1eq = +1 / (omega_RA * np.imag(Z1eq))

    # Check output quality

    if R1eq < 0. or C1eq < 0.:
        logging.warning(
            f"WARNING impedenceParallel function, There's something wrong with the calculation, R or C is negative.. R {R1eq}, C {C1eq}")

    return R1eq, C1eq


def tri2star(T1, T2, T3):
    '''Transforms three resistances in triangular connection into
    three resistances in star connection

    Parameters
    ----------
    T1 : float
        Resistance 1
    T2: float
        Resistance 2
    T3: float
        Resistance 3

    Returns
    -------
    tuple
        tuple of floats: 3 new resistances (star connection)

    Raises
    ------
    TypeError
        if not or floats
    '''

    # Check input data type

    if not isinstance(T1, float):
        raise TypeError(f'ERROR tri2star function, input T1 is not a float: T1 {T1}')
    if not isinstance(T2, float):
        raise TypeError(f'ERROR tri2star function, input T2 is not a float: T2 {T2}')
    if not isinstance(T3, float):
        raise TypeError(f'ERROR tri2star function, input T3 is not a float: T3 {T3}')

    # Check input data quality

    if T1 < 0. or T2 < 0. or T3 < 0.:
        logging.warning(
            f"WARNING tri2star funtion, There's something wrong with the input, one of them is negative.. T1 {T1}, T2 {T2}, T3 {T3}")

    T_sum = T1 + T2 + T3
    S1 = T2 * T3 / T_sum
    S2 = T1 * T3 / T_sum
    S3 = T2 * T1 / T_sum
    return S1, S2, S3


def long_wave_radiation(theta_a, SSW=1.):
    '''Estimation of sky and ground temperatures via VDI6007 model:
    theta_a outdoor air temperature [°C]
    SSW factor to count the clear non-clear sky

    Parameters
    ----------
    theta_a : numpy.array
        external temperature [°C]
    SSW : float
        factor to count the clear non-clear sky range 0-1 (1 clear sky)

    Returns
    -------
    tuple
        tupleof np.array:
        irradiance from sky vault,
        irradiance from ground,
        ground equivalent temeprature,
        sky equivalent temeperature

    Raises
    ------
    TypeError
        if not numpy array or floats
    '''

    # Check input data type

    if not isinstance(theta_a, np.ndarray):
        raise TypeError(f'ERROR longWaveRadiation function, input theta_a is not a np.array: theta_a {theta_a}')
    if not isinstance(SSW, float):
        try:
            SSW = float(SSW)
        except ValueError:
            raise TypeError(f'ERROR longWaveRadiation function, input SSW is not a float: SSW {SSW}')

            # Check input data quality

    if not np.all(np.greater(theta_a, -50.)) or not np.all(np.less(theta_a, 60.)):
        logging.warning(
            f"WARNING longWaveRadiation function, the theta_a input is outside range [-50,60]: theta_a {theta_a}")
    if not 0. <= SSW <= 1.:
        logging.warning(f"WARNING longWaveRadiation function, the SSW input is outside range [0,1]: SSW {SSW}")

    Ea_1 = 9.9 * 5.671 * 10 ** (-14) * (273.15 + theta_a) ** 6

    alpha_L = 2.30 - 7.37 * 10 ** (-3) * (273.15 + theta_a)
    alpha_M = 2.48 - 8.23 * 10 ** (-3) * (273.15 + theta_a)
    alpha_H = 2.89 - 1.00 * 10 ** (-2) * (273.15 + theta_a)

    Ea = Ea_1 * (1 + (alpha_L + (1 - (1 - SSW) / 3) * alpha_M + ((1 - (1 - SSW) / 3) ** 2) * alpha_H) * (
            (1 - SSW) / 3) ** 2.5)
    Ee = -(0.93 * 5.671 * 10 ** (-8) * (273.15 + theta_a) ** 4 + (1 - 0.93) * Ea)

    theta_erd = ((-Ee / (0.93 * 5.67)) ** 0.25) * 100 - 273.15  # [°C]
    theta_atm = ((Ea / (0.93 * 5.67)) ** 0.25) * 100 - 273.15  # [°C]

    return Ea, Ee, theta_erd, theta_atm


def loadHK(perc_rad, perc_rad_aw, perc_altro_irr, A_aw, A_raum):
    '''loadHK  -  Distribution of heat load on the nodes (surface nodes and air node) based on heat emitters
    of the building
        perc_rad: fraction of radiant heating/cooling surfaces on total heating/cooling load
        perc_rad_aw: fraction of radiant heating/cooling inside external walls on the total radiant heating/cooling load
        perc_altro_irr: this is the radiant fraction the heating/cooling systems that are not embedded in the surfaces

        Let's say that a roof has 3 systems:
            1) a radiant internal ceiling (20% of the total load)
            2) a radiant floor (30% of the total load) in contact to the ground
            3) a radiator working 80% convective and 20% radiant (50% of the total load)

        sigma_fhk: 0.5
            sum of 20% for the radiant ceiling and 30% for the radiant floor
        sigma_fhk_aw: 0.6
            because 0.3/(0.2+0.3) is equal to 0.6, i.e. the part of the radiant surface load
            associated to external surfaces
        sigma_hk_str: 0.2
            because the radiator is radiative for the 20%

    Parameters
    ----------
    perc_rad : float
        percentage of heat flow by radiant floors (sigma_fhk)
    perc_rad_aw : float
        percentage of radiant floors installed on AW walls (sigma_fhk_aw)
    perc_altro_irr : float
        percentage of radiant load (out of total which is rad+conv) by other types of heat emitters (e.g.: fan-coils, radiators) (sigma_rad_str)
    A_aw : float
        sum of the exterior opaque building components
    A_raum : float
        sum of internal partitions (all IW components) and exterior opaque building components

    Returns
    ----------
    tuple
        tuple of floats
        sigma_hk_iw: radiant heat load  on IW building components (on surface node IW)
        sigma_hk_aw: radiant heat load  on AW building components (on surface node AW)
        sigma_hk_kon: convective heat load (on air node)

    Raises
    ------
    ValueError
        if not floats
    '''

    # Check input data type

    if not isinstance(perc_rad, float):
        try:
            perc_rad = float(perc_rad)
        except ValueError:
            raise TypeError(f'ERROR loadHK function, input perc_rad is not a float: perc_rad {perc_rad}')
    if not isinstance(perc_rad_aw, float):
        try:
            perc_rad_aw = float(perc_rad_aw)
        except ValueError:
            raise TypeError(f'ERROR loadHK function, input perc_rad_aw is not a float: perc_rad_aw {perc_rad_aw}')
    if not isinstance(perc_altro_irr, float):
        try:
            perc_altro_irr = float(perc_altro_irr)
        except ValueError:
            raise TypeError(
                f'ERROR loadHK function, input perc_altro_irr is not a float: perc_altro_irr {perc_altro_irr}')
    if not isinstance(A_raum, float):
        try:
            A_raum = float(A_raum)
        except ValueError:
            raise TypeError(f'ERROR loadHK function, input A_raum is not a float: A_raum {A_raum}')
    if not isinstance(A_aw, float):
        try:
            A_aw = float(A_aw)
        except ValueError:
            raise TypeError(f'ERROR loadHK function, input A_aw is not a float: A_aw {A_aw}')

            # Check input data quality

    for p in [perc_rad, perc_rad_aw, perc_altro_irr]:
        if not 0. <= p <= 1.:
            logging.warning(
                f"WARNING loadHK function, one of the percentage inputs is outside range [0,1]: perc_rad {perc_rad}, perc_rad_aw {perc_rad_aw}, perc_altro_irr {perc_altro_irr}")
    for area in [A_raum, A_aw]:
        if not 0. <= area:
            logging.warning(
                f"WARNING loadHK function, one of the area inputs is negative: A_raum {A_raum}, A_aw {A_aw}")

            # %Note: the sum of the 3 outputs must be equal to 1

    if perc_rad == 1:

        perc_rad_iw = 1 - perc_rad_aw
        sigma_hk_iw = perc_rad_iw
        sigma_hk_aw = perc_rad_aw
        sigma_hk_kon = 0

    elif perc_rad == 0:

        sigma_hk_iw = 0
        sigma_hk_aw = 0
        sigma_hk_kon = 1

    else:
        perc_rad_iw = 1 - perc_rad_aw
        perc_altro = 1 - perc_rad
        perc_altro_irr = perc_altro * perc_altro_irr
        sigma_hk_iw = perc_altro_irr * (A_raum - A_aw) / A_raum + perc_rad_iw * perc_rad
        sigma_hk_aw = perc_altro_irr * A_aw / A_raum + perc_rad_aw * perc_rad
        sigma_hk_kon = 1 - sigma_hk_iw - sigma_hk_aw

    return [sigma_hk_iw, sigma_hk_aw, sigma_hk_kon]
