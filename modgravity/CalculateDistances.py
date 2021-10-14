import numpy as np
from scipy.integrate import quad


class CalculateDistances:
    c = 299792.458  # km/s
    h0 = 72  # km/s/Mpc
    omega_m = .3
    omega_lambda = .7

    def __init__(self):
        pass

    @classmethod
    def lum_dist(cls, z_max, z):
        """ integrates over redshift range (given by a vector z) to calculate standard luminosity distance """
        integrand = lambda z_prime: 1 / np.sqrt(cls.omega_m * (1 + z_prime) ** 3 + cls.omega_lambda)
        integral = quad(integrand, z[0], z_max)
        constant_term = cls.c * (1 + z_max) / cls.h0
        return constant_term * integral[0]

    @classmethod
    def lum_dist_array(cls, z_max, z):
        """ generates np array of luminosity distances over a uniform interval of redshifts provided by vector z """
        arr = np.array([cls.lum_dist(z_max, np.linspace(0, i)) for i in np.nditer(z)])
        for i in np.nditer(z):
            cls.lum_dist(i, np.linspace(0, i))
        return arr

    @classmethod
    def alpha_dist(cls, z_max, z, alpha):
        """ integrates over redshift range (given by a vector z) to calculate luminosity distance modified by some
        alpha """
        integrand = lambda z_prime: ((1 + z_prime) ** (alpha - 2)) / np.sqrt(cls.omega_m * (1 + z_prime) ** 3
                                                                             + cls.omega_lambda)
        integral = quad(integrand, z[0], z_max)
        constant_term = cls.c * (1 + z_max) / cls.h0
        return constant_term * integral[0]

    @classmethod
    def alpha_dist_array(cls, z_max, z, alpha):
        """ generates np array of modified alpha distances over a uniform interval of redshifts provided by vector z """
        arr = np.array([cls.lum_dist(z_max, np.linspace(0, i)) for i in np.nditer(z)])
        for i in np.nditer(z):
            cls.alpha_dist(i, np.linspace(0, i), alpha)
        return arr

    @classmethod
    def chi(cls, z_max, alpha, A_term, m_g, E_e, t_e, t_a):
        a_t_e = 1 / (1 + z_max)  # z_max corresponds to the redshift at t_e
        # a(t_e) = a(t) evaluated at t = t_e
        # given z := a0 (=1) / a(t_e) - 1 = 1 / a(t_e) - 1 => 1 / (z + 1) = a(t_e)

        def term_I(z_max, t_e, t_a):
            # a(t) = 1 / (1 + z) => 1 / a(t) = 1 + z
            return quad(lambda t: 1 + t * z_max, t_e, t_a)[0]

        def term_II(m_g, t_e, t_a):
            constant = 1 / 2 * m_g ** 2 / (a_t_e ** 2 * E_e ** 2)  # debugger value = 0.124999 (checked on wolfram)
            integral = quad(lambda t: 1 / (1 + (t * z_max)), t_e, t_a)
            return constant * integral[0]  # debugger value = 4.882812499999999e-40

        def term_III(alpha, A_term, E_e):
            constant = 1 / 2 * A_term * (a_t_e * E_e) ** (alpha - 2)
            integral = quad(lambda t: (1 / (1 + (t * z_max))) ** (1 - alpha), t_e, t_a)
            return constant * integral[0]

        return term_I(z_max, t_e, t_a) - term_II(t_e, t_a, m_g) - term_III(A_term, E_e, alpha)

    @classmethod
    def chi_array(cls, z, alpha, A_term, m_g, E_e, t_e, t_a):
        step = z[len(z) - 1] / len(z)
        z_values = np.arange(0, z[len(z) - 1], step)
        arr = np.array([cls.chi(i, m_g, E_e, t_e, t_a, alpha, A_term) for i in z_values])
        return arr
