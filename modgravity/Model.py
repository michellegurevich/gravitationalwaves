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
    def get_waveform(cls):
        return cls.pycbc_wf

    @classmethod
    def get_params(cls, wf):
        df = wf['delta_f']
        f = np.linspace(30, 200, int(170 / df))
        z_max = 1
        z = np.linspace(0, z_max)
        return f, z_max, z

    @classmethod
    def decompose_waveform(cls, wf, f, z_max, z, mass_1, mass_2):
        # calculate inner product
        chirp_mass = cls.MP.chirp_mass(wf['mass1'], wf['mass2'])
        phi_r, phi_i = cls.MP.std_polarization_array(f, z_max, z, chirp_mass, mass_1, mass_2)
        ip = np.lib.scimath.sqrt(phi_i * np.conj(phi_i))  # sqrt(-r) in R -> i*sqrt(r) in C
        phase = phi_i / ip
        return ip, phase

    @classmethod
    def perform_modification(cls, wf, std_phase, f, z_max, z, alpha, A_term):
        m_1 = wf['mass1']
        m_2 = wf['mass2']
        chirp_mass = cls.MP.chirp_mass(m_1, m_2)
        d_phase = []
        mod_phase = []

        for i in range(len(f)):
            d_phase.append(cls.MP.delta_psi(f[i], z_max, z, alpha, A_term, chirp_mass))
            mod_phase.append(std_phase + d_phase)

        return mod_phase

    @classmethod
    def mod_polarization_array(cls, f, f_cut, z_max, z, alpha, A_term):
        """ assigns product of pycbc amplitude, exp( i * modified pycbc phase) to a list, returns real and imag
        components separately """

        chirp_mass = MP.chirp_mass(cls.pycbc_wf['mass1'], cls.pycbc_wf['mass2'])
        amplitude, mod_phase = cls.perform_modification(f, z_max, z, alpha, A_term)
        arr = []

        for i in range(len(f)):
            h_tilde = amplitude * np.exp(1j * mod_phase) if f[i] < f_cut else 0
            arr.append(h_tilde)

        return [arr[i].real for i in range(len(f))], [arr[j].imag * 1j for j in range(len(f))]
