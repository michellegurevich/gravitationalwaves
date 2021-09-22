import numpy as np
import math

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from ModifiedPolarization import ModifiedPolarization
from Plots import Plots


def main():
    z_max = 4
    z = np.linspace(0, z_max)

    h0 = 72  # Hubble parameter today in km/s/Mpc
    omega_m = .3
    omega_lambda = .7
    alpha = 3
    A_term = .0001

    chirp_mass = 25e30  # kg
    f = 10e-3  # Hz

    h = 1
    f_e = 10**4  # Hz
    l_g = 1.6*10**16  # m
    E_e = h * f_e
    m_g = h / l_g
    lambda_A_term = h * A_term ** (1 / (alpha - 2))  # always has units of length irrespective of alpha value

    CD = CalculateDistances()
    # lum_dist = CD.lum_dist_array(z_max, z)
    # alpha_dist = CD.alpha_dist_array(z_max, z, 3)
    # print(lum_dist, alpha_dist)

    P = Plots()
    # P.scale_factor(z)
    # P.alpha_lum_ratios(z_max, z)
    # P.chi_to_mod_chi_ratio(z, m_g, E_e)

    MP = ModifiedPolarization()
    print(MP.delta_psi(alpha, A_term, chirp_mass, z, f))

if __name__ == '__main__':
    main()
