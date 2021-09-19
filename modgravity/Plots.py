from matplotlib import pyplot as plt
import numpy as np

from CalculateDistances import CalculateDistances
from SetCosmology import SetCosmology


class Plots:

    def __init__(self):
        pass

    # plot scale factor against redshift
    @staticmethod
    def scale_factor(z):
        SC = SetCosmology()
        hp = SC.get_hp(z)
        plt.plot(z, hp)
        plt.xlabel('$z$')
        plt.ylabel('$H(t)$')
        plt.title('$H(t)$ as a function of redshift')
        return plt.show()

    # plot ratios of (modified) alpha and (standard) luminosity distances against redshift
    @staticmethod
    def alpha_lum_ratios(z_max, z, h0, omega_m, omega_lambda):
        results = []
        CD = CalculateDistances()

        alpha_values = np.arange(0, 4.5, .5)  # alpha values end at 4, per GR tests paper specs
        for a in alpha_values:
            lum_dist_array = CD.lum_dist_array(z_max, z, h0, omega_m, omega_lambda)
            alpha_dist_array = CD.alpha_dist_array(z_max, z, a, h0, omega_m, omega_lambda)
            results = np.divide(alpha_dist_array[1:], lum_dist_array[1:])  # division by zero error, skip first element
            # results = np.insert(results, 0, 0, axis=None)
            plt.plot(np.hstack(z[1:]), np.hstack(results), label=r'$\alpha$ = '+str(a))

        plt.xlabel(r'$z$')
        plt.ylabel(r'$D_{\alpha} / D_L {Mpc}$')
        plt.title('Ratio of alpha to standard luminosity distances')
        # handles, labels = plt.gca().get_legend_handles_labels()
        # plt.gca().legend(handles=handles[::-1], labels=labels[::-1])
        plt.legend(fontsize='small')
        return plt.show()

    # plot the ratio of conformal distance chi with and without modification terms
    @staticmethod
    def chi_to_mod_chi_ratio(z_max, m_g, E_e):
        CD = CalculateDistances()

        alpha = 3
        eta_dsrt = 1 * 10e-35  # parameter of order of Planck length
        chi = CD.chi(z_max, m_g, E_e, 5, 10, 0, 0)
        mod_chi = CD.chi(z_max, m_g, E_e, 5, 10, alpha, eta_dsrt)

        # plt.plot(chi / mod_chi, z)
        # plt.xlabel(r'$z$')
        # plt.ylabel(r'$\chi$ / modified $\chi$')
        # plt.title('Ratio of standard to modified conformal distance in Double Special Relativity')

        print(chi, mod_chi)

        return 0  # plt.show
