from matplotlib import pyplot as plt
import numpy as np
import cmath
import math

from CalculateDistances import CalculateDistances
from SetCosmology import SetCosmology
from ModifiedPolarization import ModifiedPolarization


class Plots:
    CD = CalculateDistances()
    SC = SetCosmology()
    MP = ModifiedPolarization()

    def __init__(self):
        pass

    @staticmethod
    def scale_factor(z):
        """ plot scale factor against redshift """
        hp = cls.SC.get_hp(z)
        plt.plot(z, hp)
        plt.xlabel('$z$')
        plt.ylabel('$H(t)$')
        plt.title('$H(t)$ as a function of redshift')
        return plt.show()

    @staticmethod
    def alpha_lum_ratios(z_max, z):
        """ plot ratios of (modified) alpha and (standard) luminosity distances against redshift """
        alpha_values = np.arange(0, 4.5, .5)  # alpha values end at 4, per GR tests paper specs

        for a in alpha_values:
            lum_dist_array = cls.CD.lum_dist_array(z_max, z)
            alpha_dist_array = cls.CD.alpha_dist_array(z_max, z, a)
            results = np.divide(alpha_dist_array[1:], lum_dist_array[1:])  # division by zero error, skip first element
            # results = np.insert(results, 0, 0, axis=None)
            plt.plot(np.hstack(z[1:]), np.hstack(results), label=r'$\alpha$ = '+str(a))

        plt.xlabel(r'$z$')
        plt.ylabel(r'$D_{\alpha} / D_L {Mpc}$')
        plt.title('Comparison of modified luminosity distances')
        # handles, labels = plt.gca().get_legend_handles_labels()
        # plt.gca().legend(handles=handles[::-1], labels=labels[::-1])
        plt.legend(fontsize='small')
        return plt.show()

    @staticmethod
    def chi_to_mod_chi_ratio(z_max, z, m_g, E_e):
        """ plot the ratio of conformal distance chi with and without modification terms """
        alpha = 3
        eta_dsrt = 1 * 10e-35  # parameter of order of Planck length
        chi = cls.CD.chi_array(z_max, 0, 0, m_g, E_e)
        mod_chi = cls.CD.chi_array(z, alpha, .00002, m_g, E_e, 2, 10)
        print(chi, mod_chi)

        plt.plot(z, chi / mod_chi)
        plt.xlabel(r'$z$')
        plt.ylabel(r'$\chi_{e | A,\alpha=0} / \chi_e$')
        plt.title('Ratio of standard to modified conformal distance for 'r'$\alpha$'' = '
                  + str(alpha) + ', A = ' + str(A))  # Double Special Relativity')
        return plt.show()

    @staticmethod
    def modified_polarization(f, f_cut, z, z_max, alpha, A_term, chirp_mass, m_1, m_2):
        """ plot the modified polarization, h~(f), in frequency space """
        htilde_real, _ = cls.MP.mod_polarization_array(f, f_cut, z, z_max, alpha, A_term, chirp_mass, m_1, m_2)
        """ f_em / f_obs = 1 + z => f (measured as defined in paper, aka f_obs) => define array of frequency values 
        spanning the expected range for LISA """
        plt.plot(f, htilde_real)
        plt.xscale('log')  # set xscale to log for plot axes
        plt.yscale('log')  # set yscale to log for plot axes
        plt.xlabel('$lg(f)$')
        plt.ylabel(r'log($\tilde{h}$)')
        plt.title('Modified polarization in frequency space')
        return plt.show()

    @staticmethod
    def standard_polarization(f, z_max, z, m_1, m_2):
        """ plot the standard polarization, h(f), in frequency space """
        h_real, h_imag = cls.MP.std_polarization_array(f, z_max, z, m_1, m_2)
        plt.plot(f, h_real)
        plt.xscale('log')  # set xscale to log for plot axes
        plt.yscale('log')  # set yscale to log for plot axes
        plt.xlabel('$lg(f)$')
        plt.ylabel('$lg(h)$')
        plt.title('Standard polarization in frequency space')
        return plt.show()
