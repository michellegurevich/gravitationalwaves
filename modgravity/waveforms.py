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


class waveforms:

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
        cm  = lambda m, n: ((m * n) ** (3/5)) / ((m + n) ** (1/5))  # chirp mass
        eta = lambda m, n: m * n / ((m + n) ** 2)  # symmetric mass ratio
        self.chirp_mass     = cm(self.m_1, self.m_2) * SOLAR_MASS
        self.eta            = eta(self.m_1, self.m_2)

    def lambda_A(self):
        return H_PLANCK * self.A ** (1 / (self.alpha - 2))

    def u(self):
        """ mass term for amplitude """
        return math.pi * self.chirp_mass * self.f

    def beta(self):
        D_0 = self.dist.mod_luminosity(z=self.z, alpha=0)
        return (math.pi ** 2 * D_0 * self.chirp_mass) / (LAMBDA_G ** 2 * (1 + self.z[-1]))

    def zeta(self):
        if self.alpha == 1:
            D_1 = self.dist.mod_luminosity(z=self.z, alpha=1)
            zeta = (math.pi * D_1) / self.lambda_A()
        else:
            D_alpha = self.dist.mod_luminosity(z=self.z, alpha=self.alpha)
            cm = self.chirp_mass
            term_i = np.power(math.pi,(2 - self.alpha) / (1 - self.alpha))
            term_ii = D_alpha / np.power(self.lambda_A(), (2 - self.alpha))
            term_iii = (cm ** (1 - self.alpha)) / ((1 + self.z[-1]) ** (1 - self.alpha))
            zeta = term_i * term_ii * term_iii
        return zeta

    def chirp_mass_e(self, z):
        """ the (emitted) source chirp mass M_e """
        return self.chirp_mass / (1 + z[-1])

    def v_slo(self):
        """ proportional to velocity at last stable orbit in the strongly
        relativistic regime (here we use the Schwarzschild metric termination
        point) """
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

    def gr_phase(self):
        frequency_term = 2 * math.pi * F_E * T_C
        mass_term = 3 / 128 * ((self.u() ** -5 /3) * self.gr_phase_numcoeff())
        return frequency_term - PHI_C - math.pi / 4 + mass_term

    def gr_phase_numcoeff(self):
        v = self.u() ** (1/3)  # adjust mass term exponent
        gamma = .577216  # euler mascheroni constant (dimensionless)

        v_2 = 20 / 9 * (743 / 336 + (11 / 4 * self.eta))
        v_3 = 16 * math.pi
        v_4 = 10 * (3058673 / 1016064 + (5329 / 1008 * self.eta) + (617 / 144 * self.eta**2))
        v_5 = math.pi * (38645 / 756 - (65 / 9 * self.eta)) * (1 + 3 * np.log(v / self.v_slo()))
        v_6 = (11583231236531 / 4694215680
                - (640 / 3 * math.pi**2)
                - (6848 * gamma / 21)
                - (6848 / 21 * np.log(4 * v))
                + (-15737765635 / 3048192 + 2255 * (math.pi ** 2) / 12) * self.eta
                + (76055/1728 * self.eta**2)
                - (127825/1296 * self.eta**3))
        v_7 = math.pi * (77096675 / 254016 + 378515 / 1512 * self.eta - 74045 / 756 * self.eta**2)
        return 1 + (v_2 * v**2) - (v_3 * v**3) + (v_4 * v**4) + (v_5 * v**5) + (v_6 * v**6) + (v_7 * v**7)

    def delta_phase(self):
        term_i = self.beta() / self.u()
        if self.alpha == 1:
            term_ii = self.zeta() * math.log(self.u())
        else:
            term_ii = self.zeta() * (self.u() ** (self.alpha - 1))
        return -1 * (term_i + term_ii)

    def modified_phase(self):
        return self.gr_phase() + self.delta_phase()

    def h_standard(self):
        """ calculates standard wave polarization from amplitude and phase terms """
        h = self.gr_amplitude() * np.exp(1j * self.gr_phase())
        return h.real, h.imag

    def h_modified(self):
        """ calculates modified wave polarization from amplitude and phase terms """
        h_tilde = self.modified_amplitude() * np.exp(1j * self.modified_phase())
        # remove inf term from start of array (to prevent divide by 0 error)
        np.nan_to_num(h_tilde, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        return h_tilde.real, h_tilde.imag

    def decompose(self):
        # calculate inner product
        phi_r, phi_i = self.h_standard()
        ip = np.lib.scimath.sqrt(phi_i * np.conj(phi_i))  # sqrt(-r) in R -> i*sqrt(r) in C
        np.nan_to_num(ip, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        phase = phi_i / ip
        np.nan_to_num(phase, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        return ip, phase

"""

do not calculate f_dot and set equal to zero to recover max value attained by f

@classmethod
def calculate_f_max(cls, z):
    f_dot = cls.df_e_over_dt_e() / ((1 + z)**2)
    return 0

@classmethod
def df_e_over_dt_e(cls, chirp_mass, z, f_e):
    M_e = cls.chirp_mass_e(chirp_mass, z)
    m, eta = m(chirp_mass, z, m_1, m_2)
    factor = (96 / (5 * math.pi * M_e**2)) * (math.pi * M_e * f_e)**(11/3)
    # verify in debugger that eta is calc as expected
    term_ii = math.pi * m * f_e
    term_i = (743/336 + 11/4 * eta) * (term_ii**(2/3))
    return factor * (1 - term_i + 4 * math.pi * term_ii)

"""
