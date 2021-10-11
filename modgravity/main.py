import numpy as np
import math

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from ModifiedPolarization import ModifiedPolarization
from TimeDomain import TimeDomain
from Plots import Plots


def main():
    z_max = 4
    z = np.linspace(0, z_max)

    alpha = 3
    A_term = .0001
    chirp_mass = 25 * 2e30  # kg
    f = 10e-3  # Hz

    # CD = CalculateDistances()
    # lum_dist = CD.lum_dist_array(z_max, z)
    # alpha_dist = CD.alpha_dist_array(z_max, z, 3)
    # print(lum_dist, alpha_dist)

    P = Plots()
    # P.scale_factor(z)
    # P.alpha_lum_ratios(z_max, z)
    # P.chi_to_mod_chi_ratio(z, m_g, E_e)

    MP = ModifiedPolarization()
    # print(MP.delta_psi(alpha, A_term, chirp_mass, z, f))
    # print(MP.mod_polarization(f, 10e-4, chirp_mass, z, alpha, A_term))  # f_max < f WORKS FINE (returns 0)
    # print(MP.mod_polarization_log(f, 10e-2, chirp_mass, z, alpha, A_term))  # f_max > f
    # print(MP.mod_polarization_array(f, 100, chirp_mass, z, alpha, A_term))
    # P.modified_polarization(z, f, 100, chirp_mass, alpha, A_term)
    P.standard_polarization(z, f, 100, chirp_mass, 16 * 2e30, 15 * 2e30)

    TD = TimeDomain
    TD.plot_phenomA()

if __name__ == '__main__':
    main()
