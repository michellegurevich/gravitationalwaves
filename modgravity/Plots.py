from matplotlib import pyplot as plt
import numpy as np
import cmath
import math

from CalculateDistances import CalculateDistances
from SetCosmology import SetCosmology
from ModifiedPolarization import ModifiedPolarization


class Plots:

    def __init__(self):
        pass

    @staticmethod
    def scale_factor(z):
        """ plot scale factor against redshift """
        SC = SetCosmology()
        hp = SC.get_hp(z)
        plt.plot(z, hp)
        plt.xlabel('$z$')
        plt.ylabel('$H(t)$')
        plt.title('$H(t)$ as a function of redshift')
        return plt.show()

    @staticmethod
    def alpha_lum_ratios(z_max, z):
        """ plot ratios of (modified) alpha and (standard) luminosity distances against redshift """
        results = []
        CD = CalculateDistances()

        alpha_values = np.arange(0, 4.5, .5)  # alpha values end at 4, per GR tests paper specs
        for a in alpha_values:
            lum_dist_array = CD.lum_dist_array(z_max, z)
            alpha_dist_array = CD.alpha_dist_array(z_max, z, a)
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
    def chi_to_mod_chi_ratio(z, m_g, E_e):
        """ plot the ratio of conformal distance chi with and without modification terms """
        CD = CalculateDistances()

        alpha = 3
        A = .002
        eta_dsrt = 1 * 10e-35  # parameter of order of Planck length
        chi = CD.chi_array(z, m_g, E_e, 2, 10, 0, 0)
        mod_chi = CD.chi_array(z, m_g, E_e, 2, 10, alpha, .00002)
        print(chi, mod_chi)

        plt.plot(z, chi / mod_chi)
        plt.xlabel(r'$z$')
        plt.ylabel(r'$\chi_{e | A,\alpha=0} / \chi_e$')
        plt.title('Ratio of standard to modified conformal distance for 'r'$\alpha$'' = '
                  + str(alpha) + ', A = ' + str(A))  # Double Special Relativity')
        return plt.show()

    @staticmethod
    def modified_polarization(z, f, f_max, chirp_mass, alpha, A_term, m_1, m_2):
        """ plot the modified polarization, h~(f), in frequency space against redshift """
        MP = ModifiedPolarization()
        h_tilde_real, h_tilde_imag = MP.mod_polarization_array(f, f_max, chirp_mass, z, alpha, A_term, m_1, m_2)
        """ f_em / f_obs = 1 + z => f (measured as defined in paper, aka f_obs) => define array of frequency values 
        spanning the expected range for LISA """
        f_obs = np.linspace(10e-5, 10e-1)
        # plt.plot(np.log(f_obs), np.cos(h_tilde_real))
        # plt.plot(np.log(f_obs), np.sin(h_tilde_imag))
        plt.plot(np.log(f_obs), np.cos(h_tilde_real) + np.sin(h_tilde_imag))
        plt.xlabel('$lg(f)$')
        plt.ylabel(r'log($\tilde{h}$)')
        plt.title('Modified polarization in frequency space')
        return plt.show()

    @staticmethod
    def standard_polarization(z, f, f_max, chirp_mass, m_1, m_2):
        """ plot the standard polarization, h(f), in frequency and time """
        MP = ModifiedPolarization()
        h_real, h_imag = MP.std_polarization_array(f, f_max, chirp_mass, z, m_1, m_2)
        f_obs = np.linspace(10e-5, 10e-1)
        plt.plot(np.log(f_obs), np.cos(h_real) + np.sin(h_imag))
        plt.xlabel('$lg(f)$')
        plt.ylabel('$lg(h)$')
        plt.title('Standard polarization in frequency space')
        return plt.show()
