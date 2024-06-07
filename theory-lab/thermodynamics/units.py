import constants as C
from numpy import pi 

def mev_to_omega(val):
    """
    VASP output [meV] to omega [s^-1]
    """
    return val * 1e-3 * C.ELEMENTARY_CHARGE / C.RED_PLANCK_CONSTANT

def omega_to_nu(val):
    """
    Omega [s^-1] to nu [cm^-1]
    """
    return val / C.SPEED_OF_LIGHT / 2 / pi / 100
