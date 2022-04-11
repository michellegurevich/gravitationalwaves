import numpy as np
from scipy.integrate import quad


C = 299792.458  # km/s
H_PLANCK = 1
H0 = 72  # km/s/Mpc
OMEGA_M = .3
OMEGA_LAMBDA = .7

class distances:

    def __init__(self, cosmo_params, phenom_params):
        self.cosmo_params   = cosmo_params
        self.phenom_params  = phenom_params

        # unpack cosmo parameters
        self.z              = cosmo_params.get('redshift', np.linspace(0,0))
        self.f              = cosmo_params.get('frequency', np.linspace(0,0))
        self.E_e            = H_PLANCK * cosmo_params.get('f_e', 1)
        self.t_e            = cosmo_params.get('t_e', 1)
        self.t_a            = cosmo_params.get('t_a', 1)

        # unpack phenomenological parameters
        self.A              = phenom_params.get('A', 0)
        self.alpha          = phenom_params.get('alpha', 0)
        self.m_g            = H_PLANCK / phenom_params.get('lambda_g', 1)

    def __calc_luminosity(self, z):
        integrand = lambda z_prime: 1 / np.sqrt(OMEGA_M * (1 + z_prime) ** 3 + OMEGA_LAMBDA)
        integral = quad(integrand, z[0], z[-1])
        coefficient = C * (1 + z[-1]) / H0
        return coefficient * integral[0]

    def luminosity(self, z):
        return np.array([self.__calc_luminosity(np.linspace(0,z_i)) for z_i in np.nditer(z)])

    def __calc_mod_luminosity(self, z, alpha):
        integrand = lambda z_prime: ((1 + z_prime) ** (alpha - 2)) / np.sqrt(OMEGA_M * (1 + z_prime) ** 3 + OMEGA_LAMBDA)
        integral = quad(integrand, z[0], z[-1])
        coefficient = C * (1 + z[-1]) / H0
        return coefficient * integral[0]

    def mod_luminosity(self, z, alpha):
        return np.array([self.__calc_mod_luminosity(np.linspace(0,z_i), alpha) for z_i in np.nditer(z)])

    def __calc_chi_term(self, z, E_e, t_e, t_a, alpha, A, m_g):
        # a(t_e) = a(t) evaluated at t = t_e, for redshift at time of emission
        a_t_e = 1 / (1 + z)
        # z := a0 (=1) / a(t_e) - 1 = 1 / a(t_e) - 1 => 1 / (z + 1) = a(t_e)
        # a(t) = 1 / (1 + z) => 1 / a(t) = 1 + z
        term_I = quad(lambda t: 1 + t * z, t_e, t_a)[0]
        # debugger value = 0.124999 (verified on wolfram)
        II_constant = 1 / 2 * m_g ** 2 / (a_t_e ** 2 * E_e ** 2)
        II_integral = quad(lambda t: 1 / (1 + (t * z)), t_e, t_a)
        # debugger value = 4.882812499999999e-40
        term_II = II_constant * II_integral[0]
        III_constant = 1 / 2 * A * (a_t_e * E_e) ** (alpha - 2)
        III_integral = quad(lambda t: (1 / (1 + (t * z))) ** (1 - alpha), t_e, t_a)
        term_III = III_constant * III_integral[0]
        return term_I - term_II - term_III

    def chi_term(self, z):
        return np.array([self.__calc_chi_term(z_i, self.E_e, self.t_e, self.t_a, alpha=2, A=0, m_g=0) for z_i in np.nditer(z)])

    def mod_chi_term(self, z):
        return np.array([self.__calc_chi_term(z_i, self.E_e, self.t_e, self.t_a, self.alpha, self.A, self.m_g) for z_i in np.nditer(z)])
