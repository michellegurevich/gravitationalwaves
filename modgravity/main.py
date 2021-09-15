import numpy as np

from SetCosmology import SetCosmology
from Plots import Plots


def main():
    z_max = 4
    z = np.linspace(0, z_max)
    z_prime = np.linspace(0, z_max).reshape(-1, 1)
    SC = SetCosmology()
    SC.get_pars()
    SC.get_results()
    print(SC.get_bg())
    print(SC.get_ld(z))
    P = Plots()
    P.scale_factor(z)
    h0 = 72  # Hubble parameter today in km/s/Mpc
    omega_m = .3
    omega_lambda = .7
    c = 299792.458  # km /s
    P.alpha_lum_ratios(z_max)


if __name__ == '__main__':
    main()
