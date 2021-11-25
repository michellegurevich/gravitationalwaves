import camb
import numpy as np
from SetCosmology import SetCosmology
from ModifiedPolarization import ModifiedPolarization
from TimeDomain import TimeDomain


class Model:
    SC = SetCosmology()
    MP = ModifiedPolarization()
    TD = TimeDomain()

    pycbc_wf = {
        'approximant': 'TaylorF2',
        'mass1': 10, 'mass2': 15,
        'delta_f': 1.0 / 170, 'f_lower': 40,
    }

    def __init__(self):
        pass

    @classmethod
    def get_waveforms(cls):
        return 0

    @classmethod
    def decompose_waveform(cls, wf):
        # calculate inner product
        phi_r, phi_i = cls.std_polarization_array(f, z_max, np.linspace(0, z_max), m_1, m_2)
        ip = np.lib.scimath.sqrt(phi_i * np.conj(phi_i))  # sqrt(-r) in R -> i*sqrt(r) in C
        phase = phi_i / amplitude
        return ip, phase

    @classmethod
    def perform_modification(cls, f, z_max, z, alpha, A_term):
        amplitude, std_phase = cls.decompose_waveform(wf)
        mod_phase = []

        for i in range(len(f)):
            d_phase = MP.delta_psi(f[i], z_max, z, alpha, A_term)
            mod_phase.append(std_phase + d_phase[i])

        return amplitude, mod_phase

    @classmethod
    def mod_polarization_array(cls, f, f_cut, z_max, z, alpha, A_term):
        """ assigns product of pycbc amplitude, exp( i * modified pycbc phase) to a list, returns real and imag
        components separately """

        chirp_mass = MP.chirp_mass(cls.pycbc_wf['mass1'], cls.pycbc_wf['mass2'])
        amplitude, mod_phase = cls.perform_modification(f, z_max, z, alpha, A_term)
        arr = []

        for i in range(len(f)):
            h_tilde = np.log(a_tilde[i] * np.exp(1j * MP.psi(f[i], z_max, z, alpha, A_term, chirp_mass, m_1, m_2))) \
                if f[i] < f_cut else 0
            arr.append(h_tilde)

        return [arr[i].real for i in range(len(f))], [arr[j].imag * 1j for j in range(len(f))]
