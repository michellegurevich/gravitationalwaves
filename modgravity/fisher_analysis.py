import camb
import numpy as np
import matplotlib.pyplot as plt
from Model import Model


SOLAR_MASS = 4.925 * 10e-6

class fisher_analysis:

    def __init__(self, cosmo_params, phenom_params, wf_params):
        self.cosmo_params = cosmo_params
        self.phenom_params = phenom_params
        self.wf_params = wf_params

    def generate_noise(self):
        noise = []
        return noise

    def get_signal(self, *args, **kwargs):
        wf = self.wf_params
        df = self.wf_params['delta_f']
        f = np.linspace(30, 200, int(170 / df))
        z = self.cosmo_params['redshift']
        mass_1 = self.wf_params['mass1'] * SOLAR_MASS
        mass_2 = self.wf_params['mass2'] * SOLAR_MASS
        alpha = self.phenom_params['alpha']
        A = self.phenom_params['A']

        model = Model(self.cosmo_params, self.phenom_params, self.wf_params)
        amplitude, std_phase = model.decompose_waveform(wf, f, z, mass_1, mass_2)
        mod_phase = model.perform_modification(wf, std_phase, f, z[-1], z, alpha, A)
        signal = model.generate_modified_wf(wf, amplitude, mod_phase)
        return signal

    def integrate_signal(self):
        f_min = 0
        f_max = 0
        hplus = self.get_signal(wf_params['poln_cds']=[1,0])
        hcross = self.get_signal(wf_params['pln_cds']=[0,1])
        dhplus_dalpha = (hplus.amplitude() ** 2) * (e ** (1*j) * (-2) * hplus.phase()))
        return dhplus_dalpha

    def num_partial_derivs(self, cosmo_params, phenom_params, wf_params):
        # deriv calcd around some fiducial point in parameter space (will have some orienation, chirp mass, DL)
        # ex. A = alpha = 0, mg = 0
        h = self.get_signal()
        alpha = phenom_params['alpha']
        A = phenom_params['A']
        return 0

    def generate_matrix(self, cosmo_params, phenom_params, wf_params):
        # first iteration will only consider phenom param behavior

        # h+ defn of signal, ignore R term
        signal = self.get_signal()
        noise = self.generate_noise()
        Jacobian = self.num_partial_derivs()

        normalization_factor = 1 / (np.power(signal, 2) + np.power(noise, 2))
        return normalization_factor * Jacobian

    def plot(self):
        return 0 #plt.show
