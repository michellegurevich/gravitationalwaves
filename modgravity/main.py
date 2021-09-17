import numpy as np

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from Plots import Plots


def main():
    z_max = 4
    z = np.linspace(0, z_max)
    step = len(z) / z_max
    z_prime = np.linspace(0, z_max).reshape(-1, 1)
    SC = SetCosmology()
    # print(SC.get_ld(z))
    P = Plots()
    #P.scale_factor(z)

    h0 = 69  # Hubble parameter today in km/s/Mpc
    omega_m = .3
    omega_lambda = .7
    c = 299792.458  # km /s
    #P.alpha_lum_ratios(z_max, z_prime, h0, omega_m, omega_lambda)

    h = 1
    f_e = 10**4  # Hz
    l_g = 1.6*10**16  # m
    E_e = h * f_e
    m_g = h / l_g

    # try to compute integrands of lum distance from camb and own code to compare
    CD = CalculateDistances()
    print(type(SC.get_ld(z)))

    z_max = 4
    z = np.linspace(0, z_max)
    step = z_max / len(z)

    ld_results = np.array([CD.lum_dist(np.linspace(0, z), h0, omega_m, omega_lambda) for z in np.arange(0, z_max, step)])
    print(type(ld_results))

    print(SC.get_ld(z))
    print(ld_results)



    # print(CD.chi(z_max, m_g, E_e, 5, 10, 0, 0))
    # P.chi_to_mod_chi_ratio(z, z_max, m_g, E_e)


if __name__ == '__main__':
    main()
