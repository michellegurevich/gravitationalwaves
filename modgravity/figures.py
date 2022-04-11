import numpy as np
import cmath
import math
from matplotlib import pyplot as plt
from matplotlib import colors
from distances import distances
from set_cosmology import set_cosmology
from waveforms import waveforms


H_PLANCK = 1  # Planck constant, natural units

class figures:

    def __init__(self, cosmo_params, phenom_params, wf_params):
        self.cosmo_params   = cosmo_params
        self.phenom_params  = phenom_params
        self.wf_params      = wf_params
        self.dist           = distances(cosmo_params, phenom_params)
        self.sc             = set_cosmology()
        self.wf             = waveforms(cosmo_params, phenom_params, wf_params)

        # unpack cosmo parameters
        self.z              = cosmo_params.get('redshift', np.linspace(0,0))
        self.f              = cosmo_params.get('frequency', np.linspace(0,0))
        self.E_e            = H_PLANCK * cosmo_params.get('f_e', 1)
        self.t_e            = cosmo_params.get('t_e', 1)
        self.t_a            = cosmo_params.get('t_a', 1)

        # unpack phenomenological parameters
        self.A              = phenom_params.get('A', 0)
        self.alpha          = phenom_params.get('alpha', 0)
        self.m_g            = H_PLANCK / phenom_params.get('lambda_g', 1)

        # unpack waveform parameters
        self.approximant    = wf_params.get('TaylorF2', None)
        self.m_1            = wf_params.get('mass1', None)
        self.m_2            = wf_params.get('mass2', None)
        self.df             = wf_params.get('delta_f', None)
        self.f_min          = wf_params.get('f_lower', None)

    def scale_factor(self):
        """ plot scale factor against redshift """
        hp = self.sc.get_hp(self.z)
        plt.plot(self.z, hp)
        plt.xlabel('$z$')
        plt.ylabel('$H(t)$')
        plt.title('$H(t)$ as a function of redshift')
        return plt.show()

    def alpha_distance_ratio(self):
        """ plot ratios of (modified) alpha and (standard) luminosity distances against redshift """
        alpha_values = np.arange(0, 4.5, .5)  # alpha values end at 4, per GR tests paper specs
        for a in alpha_values:
            D_l = self.dist.luminosity(self.z)
            D_a = self.dist.mod_luminosity(self.z, a)
            results = np.divide(D_l[1:], D_a[1:])  # division by zero error, skip first element
            # results = np.insert(results, 0, 0, axis=None)
            plt.plot(np.hstack(np.linspace(0,12)[1:]), np.hstack(results), label=r'$\alpha$ = '+str(a))#, c=colors.Colormap('coolwarm'))

        plt.xlabel(r'$z$')
        plt.ylabel(r'$D_L / D_{\alpha}$ Mpc')
        plt.title('Modified luminosity distance')
            # handles, labels = plt.gca().get_legend_handles_labels()
            # plt.gca().legend(handles=handles[::-1], labels=labels[::-1])
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        return plt.show()

    def chi_ratio(self):
        """ plot the ratio of conformal distance chi with and without modification terms """
        eta_dsrt = 1 * 10e-35  # parameter of order Planck length
        chi = self.dist.chi_term(self.z)
        mod_chi = self.dist.mod_chi_term(self.z)

        plt.plot(self.z, chi / mod_chi)
        plt.xlabel(r'$z$')
        plt.ylabel(r'$\chi_{e | A,\alpha=0} / \chi_e$')
        plt.title('Ratio of standard to modified conformal distance for 'r'$\alpha$'' = '
                  + str(self.alpha) + ', A = ' + str(self.A))  # Double Special Relativity')
        return plt.show()

    def h_standard(self):
        """ plot the standard polarization, h(f), in frequency space """
        h_plus, h_cross = self.wf.h_standard()
        plt.plot(self.f, h_plus)
        #plt.xscale('log')  # set xscale to log for plot axes
        #plt.yscale('log')  # set yscale to log for plot axes
        plt.xlabel('$log(f)$')
        plt.ylabel('$log(h)$')
        plt.title('Standard polarization in frequency space')
        return plt.show()

    def h_modified(self):
        """ plot the modified polarization, h~(f), in frequency space """
        h_plus, _ = self.wf.h_modified()
        """ f_em / f_obs = 1 + z => f (measured as defined in paper, aka f_obs) => define array of frequency values
        spanning the expected range for LISA """
        plt.plot(self.f, h_plus)
        #plt.xscale('log')  # set xscale to log for plot axes
        #plt.yscale('log')  # set yscale to log for plot axes
        plt.xlabel('$log(f)$')
        plt.ylabel(r'log($\tilde{h}$)')
        plt.title('Modified polarization in frequency space')
        return plt.show()

    def phase_check(self, phi_r, phi_i):
        # calculate inner product
        ip, phase = self.wf.decompose()

        # plot A(f) - which should equal inner product
        plt.subplot(2, 1, 1)
        print(len(self.f), len(ip))
        plt.plot(self.f, ip)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')

        # plot e^(i*Psi) - calculated as phi / mag(phi)
        plt.subplot(2, 1, 2)
        plt.plot(self.f, phase)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Phase')
        return plt.show()
