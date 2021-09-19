import numpy as np
import math

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from Plots import Plots


def main():
    z_max = 4
    z = np.linspace(0, z_max)

    h0 = 72  # Hubble parameter today in km/s/Mpc
    omega_m = .3
    omega_lambda = .7
    alpha = 3

    h = 1
    f_e = 10**4  # Hz
    l_g = 1.6*10**16  # m
    E_e = h * f_e
    m_g = h / l_g

    CD = CalculateDistances()
    lum_dist = CD.lum_dist_array(z_max, z, h0, omega_m, omega_lambda)
    alpha_dist = CD.alpha_dist_array(z_max, z, alpha, h0, omega_m, omega_lambda)
    print(lum_dist, alpha_dist)

    chi_e = CD.chi(z_max, m_g, E_e, 2, 10, 3, 10e-35)
    print(chi_e)

    P = Plots()
    P.scale_factor(z)

if __name__ == '__main__':
    main()
