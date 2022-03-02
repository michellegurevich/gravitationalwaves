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

    def generated_noise(self):
        noise = []
        return noise

    def integrated_signal(self):
        self.wf_params.set_poln_cds((wf_params['poln_cds']=[1,0]))
        self.wf_params.set_poln_cds((wf_params['poln_cds']=[0,1]))
        h_plus, h_cross = self.wf.h_modified()
        integrand = lambda f_prime: 1 + 5
        integral = quad(integrand, self.f[0], self.f[-1])
        return

    def partial_derivs(self):
        """ deriv calculatedd around some fiducial point in parameter space (will
        have some orienation, chirp mass, distances) ex. A = alpha = m_g = 0 """
        h = self.get_signal()
        dwf_dalpha  = (hplus.amplitude() ** 2) * (e ** (1*j) * (-2) * hplus.phase()))
        return dhplus_dalpha

    def generate_matrix(self):
        """ first iteration will only consider phenom param behavior and h+ defn of
        signal, ignore R term """
        signal = self.integrated_signal()
        noise = self.generated_noise()
        Jacobian = self.partial_derivs()

        normalization_factor = 1 / (np.power(signal, 2) + np.power(noise, 2))
        return normalization_factor * Jacobian
