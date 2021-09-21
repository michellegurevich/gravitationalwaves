import numpy as np
from scipy.integrate import simps
from scipy.integrate import quad
from sympy import lambdify


class CalculateDistances:

    def __init__(self):
        pass

    @classmethod
    def lum_dist(cls, z, h0, omega_m, omega_lambda):
        c = 299792.458  # km /s

        integrand = lambda z_prime: 1 / np.sqrt(omega_m * (1 + z_prime) ** 3 + omega_lambda)
        integral = quad(integrand, z[0], z[len(z) - 1])
        constant_term = c * (1 + z[len(z) - 1]) / h0  # FROM eqn (2) THERE IS FACTOR OF POWERS OF C IN ENERGY TERM
        return constant_term * integral[0]

    @classmethod
    def lum_dist_array(cls, z_max, z, h0, omega_m, omega_lambda):
        step = z_max / len(z)
        z_values = np.arange(0, z_max, step)
        arr = np.array([cls.lum_dist(np.linspace(0, z), h0, omega_m, omega_lambda) for z in z_values])
        return arr

    @classmethod
    def alpha_dist(cls, z, alpha, h0, omega_m, omega_lambda):
        c = 299792.458  # km /s

        integrand = lambda z_prime: ((1 + z_prime) ** (alpha - 2)) / np.sqrt(
            omega_m * (1 + z_prime) ** 3 + omega_lambda)
        integral = quad(integrand, z[0], z[len(z) - 1])
        constant_term = c * (1 + z[len(z) - 1]) / h0  # FROM eqn (2) THERE IS FACTOR OF POWERS OF C IN ENERGY TERM
        return constant_term * integral[0]

    @classmethod
    def alpha_dist_array(cls, z_max, z, alpha, h0, omega_m, omega_lambda):
        step = z_max / len(z)
        z_values = np.arange(0, z_max, step)
        arr = np.array([cls.alpha_dist(np.linspace(0, z), alpha, h0, omega_m, omega_lambda) for z in z_values])
        return arr

    @classmethod
    def chi(cls, z_max, m_g, E_e, t_e, t_a, alpha, A_term):
        a_t_e = 1 / (1 + z_max)  # z_max corresponds to the redshift at t_e

        # a(t_e) = a(t) evaluated at t = t_e
        # given z := a0 (=1) / a(t_e) - 1 = 1 / a(t_e) - 1 => 1 / (z + 1) = a(t_e)

        def term_I(z_max, t_e, t_a):
            # a(t) = 1 / (1 + z) => 1 / a(t) = 1 + z
            return quad(lambda t: 1 + t * z_max, t_e, t_a)[0]

        def term_II(t_e, t_a, m_g):
            constant = 1 / 2 * m_g ** 2 / (a_t_e ** 2 * E_e ** 2)  # debugger value = 0.124999 (checked on wolfram)
            integral = quad(lambda t: 1 / (1 + (t * z_max)), t_e, t_a)
            return constant * integral[0]  # debugger value = 4.882812499999999e-40

        def term_III(A_term, E_e, alpha):
            constant = 1 / 2 * A_term * (a_t_e * E_e) ** (alpha - 2)
            integral = quad(lambda t: (1 / (1 + (t * z_max))) ** (1 - alpha), t_e, t_a)
            return constant * integral[0]

        return term_I(z_max, t_e, t_a) - term_II(t_e, t_a, m_g) - term_III(A_term, E_e, alpha)

    @classmethod
    def chi_array(cls, z, m_g, E_e, t_e, t_a, alpha, A_term):
        step = z[len(z) - 1] / len(z)
        z_values = np.arange(0, z[len(z) - 1], step)
        arr = np.array([cls.chi(i, m_g, E_e, t_e, t_a, alpha, A_term) for i in z_values])
        return arr
