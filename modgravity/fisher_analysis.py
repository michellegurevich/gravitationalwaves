import numpy as np
from scipy.integrate import quad
from scipy.fft import fft
from waveforms import waveforms

from waveforms import waveforms


SOLAR_MASS = 4.925 * 10e-6

class fisher_analysis:

    def __init__(self, cosmo_params, phenom_params, wf_params):
        self.cosmo_params   = cosmo_params
        self.phenom_params  = phenom_params
        self.wf_params      = wf_params
        self.wf             = waveforms(cosmo_params, phenom_params, wf_params)
        self.f_min          = 10e0 # Hz
        self.f_max          = 10e4 # Hz

    def generate_noise(self):
        f = self.cosmo_params['frequency']
        # comptuted for Adv. LIGO like detector
        f_0 = 215 # Hz
        f_s = 20 # Hz
        S_0 = 10e-49 # 1/Hz
        x = f / f_0
        spectral_density = np.zeros(shape=(1, len(f)))

        if f > f_s:
            for i in range(f):
                term_i   = 10 ** (16-4*(f[i]-7.9)**2)
                term_ii  = 2.4 * 10**(-62) * x**(-50)
                term_iii = 0.08 * x**(-4.69)
                term_iv  = 123.35 * (( 1 - 0.23*x**2 + 0.0764 * x**4 ) / ( 1 + 0.17 * x**2 ))
                spectral_density[i] = S_0(term_i + term_ii + term_iii + term_iv)
        return spectral_density, x

    def integrated_signal(self):
        # get noise spectral density
        sd, _ = self.generate_noise()

        # retrieve modified waveform and combine real and imaginary values element-wise
        h_tilde_re, h_tilde_im = self.wf.h_modified()
        h_tilde = h_tilde_re + h_tilde_im

        # calculate complex conjugate and associated inner product element-wise
        integrand = lambda f_prime: 2 / sd * (np.conjugate(h_tilde) * h_tilde)
        integral = quad(integrand, self.f_min, self.f_max)
        hh = 2 * integral
        snr = np.sqrt(hh)
        return snr

    def partial_derivs(self):
        """ deriv calculatedd around some fiducial point in parameter space (will
        have some orienation, chirp mass, distances) eg. zeta = beta = 0 """
        _, x = self.generate_noise()
        h_tilde_re, h_tilde_im = self.wf.h_modified()
        h_tilde = h_tilde_re + h_tilde_im

        d_ln_amp_s = h_tilde
        d_psi_s = -1j * h_tilde
        d_beta = -1j * (1 / self.wf.u()) * h_tilde
        d_t_c = 2 * math.pi * 1j * x * h_tilde

        partials = np.array(d_ln_amp_s, d_psi_s, d_beta, d_t_c)
        return partials

        # the more frequencies you have the more steep the gradients for the fisher
        # derivatives; number of events Nx initially use a constant and say that all the
        # parameters are the same; leading coefficient in front of the integrated signal

    def generate_matrix(self):
        noise, _ = self.generated_noise()
        signal = self.integrated_signal()
        partials = self.partial_derivs()

        normalization_factor = 1 / (np.power(signal, 2) + np.power(noise, 2))

        matrix = np.zeros(shape=(4, 4))
        for i, j in range(len(matrix)):
            matrix[i, j] = partials[i] * partials[j]
        return normalization_factor * matrix
