import numpy as np
from scipy.integrate import simps
from scipy.integrate import quad
from sympy import lambdify


class CalculateDistances:

    def lum_dist(self, z, h0, omega_m, omega_lambda):
        integral = quad(lambda z_prime: 1 / np.sqrt(omega_m * (1 + z_prime) ** 3 + omega_lambda), z[0], z[len(z)-1])
        c = 299792.458  # km /s
        constant_term = c * (1 + 4) / h0  # FROM eqn (2) THERE IS FACTOR OF POWERS OF C IN ENERGY TERM
        # constant_term = (1 + z_max) / h0

        return constant_term * integral[0]  # graviton travels from source to observer

    def alpha_dist(self, a, z_max, z_prime, h0, omega_m, omega_lambda):
        def integrand(z_prime, alpha):
            n = np.array(range(0, 50))  # controls number of points integrand is evaluated at
            term = omega_m * (1 + n * z_prime) ** 3 + omega_lambda
            return ((1 + z_prime) ** (alpha - 2)) / np.sqrt(term)

        integral = simps(integrand(z_prime, a), z_prime, axis=0)
        constant_term = (1 + z_max) ** (1 - a) / h0

        return constant_term * integral

    def chi(self, z_max, m_g, E_e, t_e, t_a, alpha, A_term):
        a_t_e = 1 / (1 + z_max)  # z_max corresponds to the redshift at t_e
        # a(t_e) = a(t) evaluated at t = t_e
        # given z := a0 (=1) / a(t_e) - 1 = 1 / a(t_e) - 1 => 1 / (z + 1) = a(t_e)

        def term_I(self, z_max, t_e, t_a):
            # a(t) = 1 / (1 + z) => 1 / a(t) = 1 + z
            return quad(lambda t: 1 + t * z_max, t_e, t_a)[0]

        def term_II(self, t_e, t_a, m_g):
            constant = 1/2 * m_g**2 / (a_t_e**2 * E_e**2)  # debugger value = 0.124999 (checked on wolfram)
            integral = quad(lambda t: 1 / (1 + (t * z_max)), t_e, t_a)
            return constant * integral[0]  # debugger value = 4.882812499999999e-40

        def term_III(self, A_term, E_e, alpha):
            constant = 1/2 * A_term * (a_t_e * E_e)**(alpha-2)
            integral = quad(lambda t: (1 / (1 + (t * z_max)))**(1-alpha), t_e, t_a)
            return constant * integral[0]

        return term_I(self, z_max, t_e, t_a) - term_II(self, t_e, t_a, m_g) - term_III(self, A_term, E_e, alpha)

