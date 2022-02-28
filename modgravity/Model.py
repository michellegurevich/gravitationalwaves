import camb
import numpy as np
from set_cosmology import set_cosmology
from ModifiedPolarization import ModifiedPolarization


class Model:
    bg = set_cosmology()
    

    def __init__(self, cosmo_params, phenom_params, wf_params):
        self.cosmo_params = cosmo_params
        self.phenom_params = phenom_params
        self.wf_params = wf_params
        self.MP = ModifiedPolarization()

    def decompose_waveform(self, wf, f, z, mass_1, mass_2):
        # calculate inner product
        chirp_mass = self.MP.chirp_mass(wf['mass1'], wf['mass2'])
        phi_r, phi_i = self.MP.std_polarization_array(f, z[-1], z, chirp_mass, mass_1, mass_2)
        ip = np.lib.scimath.sqrt(phi_i * np.conj(phi_i))  # sqrt(-r) in R -> i*sqrt(r) in C
        phase = phi_i / ip
        return ip, phase

    def perform_modification(self, wf, std_phase, f, z_max, z, alpha, A_term):
        m_1 = wf['mass1']
        m_2 = wf['mass2']
        chirp_mass = self.MP.chirp_mass(m_1, m_2)
        d_phase = []
        mod_phase = []

        #for i in range(len(f)):
        #    d_phase.append(self.MP.delta_psi(f[i], z_max, z, alpha, A_term, chirp_mass))
        #    mod_phase.append(std_phase + d_phase)

        return mod_phase

    @classmethod
    def generate_modified_wf(cls, wf, amplitude, mod_phase):
        modified_wf = [] # USE MOD_POL array below to recombine
        return wf['approximant'], modified_wf

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
