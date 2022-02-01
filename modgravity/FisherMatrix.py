import camb
import numpy as np
from Model import Model


class FisherMatrix:
    M = Model()

    def __init__(self):
        pass

    @classmethod
    def generate_noise(cls):
        noise = []
        return noise

    @classmethod
    def get_signal(cls, alpha, A_term):
        ''' add pass through for other parameters for modification '''
        standard_parameters = []
        modified_parameters = []
        # pycbc may use an n set of input parameters and will want those passed in together

        wf = M.get_waveform()
        f, z_max, z = M.get_params(wf)
        amplitude, std_phase = M.decompose_waveform(wf, f, z_max, z, wf['mass1'] * 4.925 * 10e-6,
                                                    wf['mass2'] * 4.925 * 10e-6)
        mod_phase = M.perform_modification(wf, std_phase, f, z_max, z, alpha, A_term)
        signal = M.generate_modified_wf(wf, amplitude, mod_phase)
        return signal

    @classmethod
    def num_partial_derivs(cls):
        # deriv calcd around some fiducial point in parameter space (will have some orienation, chirp mass, DL)
        # ex. A = alpha = 0, mg = 0
        return 0

    @classmethod
    def generate_matrix(cls):
        signal = cls.get_signal(alpha=4, A_term=.001)
        noise = cls.generate_noise()
        Jacobian = cls.num_partial_derivs

        # assuming inputs for signal and output are arrays
        normalization_factor = 1 / (np.power(signal, 2) + np.power(noise, 2))
        return normalization_factor * Jacobian

    @classmethod
    def plot(cls):
        return 0
