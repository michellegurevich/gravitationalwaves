import numpy as np
from scipy.integrate import simps
from scipy.integrate import quad
from sympy import lambdify


class CalculateDistances:

    def lum_dist(self, z_max, z_prime, h0, omega_m, omega_lambda):
        # FROM eqn (2) THERE IS FACTOR OF POWERS OF C IN ENERGY TERM ->

        def integrand(z_prime):
            n = np.array(range(0, 50))  # controls number of points integrand is evaluated at
            term = omega_m * (1 + n * z_prime) ** 3 + omega_lambda
            return 1 / np.sqrt(term)

        integral = simps(integrand(z_prime), z_prime, axis=0)
        # constant_term = c * (1 + z_max) / h0
        constant_term = (1 + z_max) / h0

        return constant_term * integral  # graviton travels from source to observer

    def alpha_dist(self, a, z_max, z_prime, h0, omega_m, omega_lambda):
        def integrand(z_prime, alpha):
            n = np.array(range(0, 50))  # controls number of points integrand is evaluated at
            term = omega_m * (1 + n * z_prime) ** 3 + omega_lambda
            return ((1 + z_prime) ** (alpha - 2)) / np.sqrt(term)

        integral = simps(integrand(z_prime, a), z_prime, axis=0)
        constant_term = (1 + z_max) ** (1 - a) / h0

        return constant_term * integral

    def chi(self, z_max, m_g, E_e, t_e, t_a, alpha, A_term):
        a_t = 1 / (1 + z_max)

        def term_I(self, t_e, t_a):
            # a(t) = 1 / (1 + z) => 1 / a(t) = 1 + z
            f = lambda t: 1 / a_t
            return quad(f, t_e, t_a)

        def term_II(self, a, t_e, t_a, m_g):
            return 0

        def term_III(self, A_term, a, E_e, alpha):
            return 0

        chi = term_I(self, t_e, t_a)

        return chi
