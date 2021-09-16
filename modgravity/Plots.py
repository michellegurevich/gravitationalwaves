from matplotlib import pyplot as plt
import numpy as np

from CalculateDistances import CalculateDistances
from SetCosmology import SetCosmology


class Plots:

    # plot scale factor against redshift
    def scale_factor(self, z):
        SC = SetCosmology()
        hp = SC.get_hp(z)
        plt.plot(z, hp)
        plt.xlabel('$z$')
        plt.ylabel('$H(t)$')
        plt.title('$H(t)$ as a function of redshift')
        return plt.show()

    # plot ratios of (modified) alpha and (standard) luminosity distances against redshift
    def alpha_lum_ratios(self, z_max, z_prime, h0, omega_m, omega_lambda):
        results = []
        CD = CalculateDistances()

        alpha_values = np.arange(0, 4.5, .5)  # alpha values end at 4, per GR tests paper specs
        for a in alpha_values:
            alpha_lum_ratio = (CD.alpha_dist(a, z_max, z_prime, h0, omega_m, omega_lambda)
                               / CD.lum_dist(z_max, z_prime, h0, omega_m, omega_lambda))[1]
            results.append(alpha_lum_ratio)

        plt.plot(np.hstack(alpha_values), np.hstack(results))
        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'$D_{\alpha} / D_L {Mpc}$')
        plt.title('Ratio of alpha to standard luminosity distances')
        return plt.show()

    # plot the ratio of conformal distance chi with and without modification terms
    def chi_to_mod_chi_ratio(self, z, z_max, m_g, E_e):
        CD = CalculateDistances()

        alpha = 3
        eta_dsrt = 1 * 10e-35  # parameter of order of Planck length
        chi = CD.chi(z_max, m_g, E_e, 5, 10, 0, 0)
        mod_chi = CD.chi(z_max, .0001, .001, 5, 10, alpha, eta_dsrt)

        #plt.plot(chi / mod_chi, z)
        #plt.xlabel(r'$z$')
        #plt.ylabel(r'$\chi$ / modified $\chi$')
        #plt.title('Ratio of standard to modified conformal distance in Double Special Relativity')

        print(chi, mod_chi)

        return 0#plt.show