import math
import cmath
import numpy as np

from matplotlib import pyplot as plt
from distances import distances


H_PLANCK = 1  # Planck constant, natural units
SOLAR_MASS = 4.925 * 10e-6
F_E = 10 ** 4  #Hz
LAMBDA_G = 1.6 * 10 ** 16  # m
EPSILON_LIGO = 1  # constant term for LIGO
# EPSILON_LISA = math.sqrt(3) / 2  # for LISA
T_C = 1  # verify in references
PHI_C = .0001  # verify in references


class modify_waveform:

    def __init__(self, cosmo_params, phenom_params, wf_params):
        self.cosmo_params   = cosmo_params
        self.phenom_params  = phenom_params
        self.wf_params      = wf_params
        self.dist           = distances(cosmo_params, phenom_params)

        # unpack cosmo parameters
        self.z      = cosmo_params.get('redshift', np.linspace(0,0))
        self.f      = cosmo_params.get('frequency', np.linspace(0,0))
        self.E_e    = H_PLANCK * cosmo_params.get('f_e', 1)
        self.t_e    = cosmo_params.get('t_e', 1)
        self.t_a    = cosmo_params.get('t_a', 1)

        # unpack phenomenological parameters
        self.A      = phenom_params.get('a', 0)
        self.alpha  = phenom_params.get('alpha', 0)
        self.m_g    = H_PLANCK / phenom_params.get('lambda_g', 1)

        # unpack waveform parameters
        self.approximant    = wf_params.get('TaylorF2', None)
        self.m_1            = wf_params.get('mass1', None)
        self.m_2            = wf_params.get('mass2', None)
        self.df             = wf_params.get('delta_f', None)
        self.f_min          = wf_params.get('f_lower', None)

        # lambda functions
        cm = lambda m, n: ((m * n) ** (3/5)) / ((m + n) ** (1/5))  # chirp mass
        eta = lambda m, n: m * n / ((m + n) ** 2)  # symmetric mass ratio
        self.chirp_mass     = cm(self.m_1, self.m_2) * SOLAR_MASS
        self.m              = eta(self.m_1, self.m_2)

    def lambda_A(self):
        return h_planck * self.A ** (1 / (self.alpha - 2))

    def u(self):
        """ mass term for amplitude """
        return math.pi * self.chirp_mass * self.f

    def beta(self):
        D_0 = self.dist.mod_luminosity(z=self.z, alpha=0)
        return (math.pi ** 2 * D_0 * self.chirp_mass) / (LAMBDA_G ** 2 * (1 + self.z[-1]))

    def zeta(self):
        lambda_A_term = self.lambda_A_term(self.alpha, self.A_term)
        if self.alpha == 1:
            D_1 = self.dist.mod_luminosity(z=self.z, alpha=1)
            zeta = (math.pi * D_1) / lambda_A_term
        else:
            D_alpha = self.dist.mod_luminosity(z=self.z, alpha=self.alpha)
            cm = self.chirp_mass
            term_i = math.pi ** (2 - self.alpha) / (1 - self.alpha)
            term_ii = D_alpha / (lambda_A_term ** (2 - self.alpha))
            term_iii = (cm ** (1 - self.alpha)) / ((1 + self.z[-1]) ** (1 - self.alpha))
            zeta = term_i * term_ii * term_iii
        return zeta

    def chirp_mass_e(self, z):
        """ the (emitted) source chirp mass M_e """
        return self.chirp_mass / (1 + z[-1])

    def v_slo(self):
        """ proportional to velocity at last stable orbit in the strongly relativistic regime (here we use the
        Schwarzschild metric termination point) """
        return 1 / math.sqrt(6)

    def gr_amplitude(self):
        """ curved A term; has frequency dependence """
        # luminosity distance calculated by skipping first redshift value z_1 = 0
        # D_L = self.dist.luminosity(z=self.z[1:])
        D_L = self.dist.luminosity(z=self.z)
        return math.sqrt(math.pi / 30) * (self.chirp_mass ** 2 / D_L)

    def modified_amplitude(self):
        """ A~(f) term """
        return EPSILON_LIGO * self.gr_amplitude() * (self.u() ** (-7 / 6))

    def gr_phase():
        return 0

    def gr_phase_numerical_coefficients():
        return 0

    def phase_perturbation():
        return 0

    def modified_phase():
        return 0

    def standard_polarization():
        return 0
