import numpy as np

from SetCosmology import SetCosmology
from Plots import Plots


def main():
    z_max = 4
    z = np.linspace(0, z_max)
    z_prime = np.linspace(0, z_max).reshape(-1, 1)
    SC = SetCosmology()
    print(SC.get_ld(z))
    P = Plots()
    # P.scale_factor(z)

    h0 = 72  # Hubble parameter today in km/s/Mpc
    omega_m = .3
    omega_lambda = .7
    c = 299792.458  # km /s
    P.alpha_lum_ratios(z_max, z_prime, h0, omega_m, omega_lambda)


if __name__ == '__main__':
    main()
