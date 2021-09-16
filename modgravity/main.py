import numpy as np

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from Plots import Plots


def main():
    z_max = 4
    # z = np.linspace(0, z_max)
    z_prime = np.linspace(0, z_max).reshape(-1, 1)
    SC = SetCosmology()
    # print(SC.get_ld(z))
    P = Plots()
    # P.scale_factor(z)

    h0 = 72  # Hubble parameter today in km/s/Mpc
    omega_m = .3
    omega_lambda = .7
    c = 299792.458  # km /s
    # P.alpha_lum_ratios(z_max, z_prime, h0, omega_m, omega_lambda)

    # try to compute integrands of lum distance from camb and own code to compare
    CD = CalculateDistances()
    # print(SC.get_ld(z) * (h0/(1+z_max)))
    # print(CD.lum_dist(z_max, z_prime, h0, omega_m, omega_lambda) * (h0/(1+z_max)))
    print(CD.chi(z_max, .0001, .001, 5, 10, 0, 0))


if __name__ == '__main__':
    main()
