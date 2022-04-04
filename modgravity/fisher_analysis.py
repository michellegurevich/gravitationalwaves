import camb
import numpy as np
from scipy.integrate import quad

from waveforms import waveforms


SOLAR_MASS = 4.925 * 10e-6

class fisher_analysis:

    def __init__(self, cosmo_params, phenom_params, wf_params):
        self.cosmo_params   = cosmo_params
        self.phenom_params  = phenom_params
        self.wf_params      = wf_params
        self.wf             = waveforms(cosmo_params, phenom_params, wf_params)

    def generate_noise(self):
        xf = fft(x - x.mean())
        Sxx = 2 * dt ** 2 / T * (xf * xf.conj())  # Compute spectrum
        Sxx = Sxx[:int(len(x) / 2)]
        return noise

    def integrated_signal(self):
        self.wf_params.set_poln_cds()
        self.wf_params.set_poln_cds()
        h_plus, h_cross = self.wf.h_modified()
        integrand = lambda f_prime: 1 + 5
        integral = quad(integrand, self.f_min, self.f_max)
        return

        # N(f) and partial derivatives ds/dalpha(f)

    def partial_derivs(self):
        """ deriv calculatedd around some fiducial point in parameter space (will
        have some orienation, chirp mass, distances) ex. A = alpha = m_g = 0 """
        dwf_dalpha  = ((hplus.amplitude() ** 2) * (e ** (1*j) * (-2) * hplus.phase()))
        return dhplus_dalpha

        # the more frequencies you have the more steep the gradients for the fisher
        # derivatives
        # number of events Nx initially use a constant and say that all the
        # parameters are the same
        # leading coefficient in front of the integrated signal

    def generate_matrix(self):
        """ first iteration will only consider phenom param behavior and h+ defn of
        signal, ignore R term """
        signal = self.integrated_signal()
        noise = self.generated_noise()
        Jacobian = self.partial_derivs()

        normalization_factor = 1 / (np.power(signal, 2) + np.power(noise, 2))
        return normalization_factor * Jacobian
