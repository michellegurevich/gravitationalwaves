import math
import cmath
import numpy as np

from CalculateDistances import CalculateDistances


class ModifiedPolarization:
    h = 1  # Planck constant
    f_e = 10 ** 4  # Hz
    lambda_g = 1.6 * 10 ** 16  # m
    epsilon = math.sqrt(3) / 2  # for LISA, for LIGO would equal 1

    E_e = h * f_e
    m_g = h / lambda_g
    t_c = 1
    phi_c = 1

    h0 = 72
    omega_m = .3
    omega_lambda = .7

    def __init__(self):
        pass

    @classmethod
    def lambda_A_term(cls, A_term, alpha):
        return cls.h * A_term ** (1 / (alpha - 2))  # always has units of length irrespective of alpha value

    @classmethod
    def u(cls, chirp_mass, f):
        return math.pi * chirp_mass * f

    @classmethod
    def beta(cls, chirp_mass, z):
        CD = CalculateDistances()
        D_0 = CD.alpha_dist(z, 0)
        return (math.pi ** 2 * D_0 * chirp_mass) / (cls.lambda_g ** 2 * (1 + z))

    @classmethod
    def zeta(cls, z, alpha, A_term, chirp_mass):
        lambda_A_term = cls.lambda_A_term(A_term, alpha)
        CD = CalculateDistances()
        if alpha == 1:
            D_1 = CD.alpha_dist(z, 1)
            zeta = (math.pi * D_1) / lambda_A_term
        else:
            D_alpha = CD.alpha_dist(z, alpha)
            term_i = math.pi ** (2 - alpha) / (1 - alpha)
            term_ii = D_alpha / (lambda_A_term ** (2 - alpha))
            term_iii = (chirp_mass ** (1 - alpha)) / ((1 + z) ** (1 - alpha))
            zeta = term_i * term_ii * term_iii
        return zeta

    @classmethod
    def delta_psi(cls, alpha, A_term, chirp_mass, z, f):
        term_i = cls.beta(chirp_mass, z) / cls.u(chirp_mass, f)
        if alpha == 1:
            term_ii = cls.zeta(z, alpha, A_term, chirp_mass) * math.ln(cls.u(chirp_mass, f))
        else:
            term_ii = cls.zeta(z, alpha, A_term, chirp_mass) * (cls.u(chirp_mass, f) ** (alpha - 1))
        return -1 * (term_i + term_ii)

    @classmethod
    def curved_A(cls, chirp_mass, z):
        CD = CalculateDistances()
        D_L = CD.lum_dist(z)
        return math.sqrt(math.pi / 30) * (chirp_mass ** 2 / D_L)

    @classmethod
    def mod_amplitude(cls, chirp_mass, z, f):
        """ tilde A(f) term """
        return cls.epsilon * cls.curved_A(chirp_mass, z) * (cls.u(chirp_mass, f) ** (-7 / 6))

    @classmethod
    def psi_gr(cls, f, chirp_mass):  # , t_c, phi_c):
        freq_term = 2 * math.pi * f * cls.t_c
        numerical_term = 3 / 128 * (cls.u(chirp_mass, f) ** -5 / 3)  # ADD CROSS PRODUCT
        return freq_term - cls.phi_c - math.pi / 4 + numerical_term

    @classmethod
    def psi(cls, alpha, A_term, chirp_mass, z, f):
        psi_gr = cls.psi_gr(f, chirp_mass)  # , t_c, phi_c)
        delta_psi = cls.delta_psi(alpha, A_term, chirp_mass, z, f)
        return psi_gr + delta_psi

    # @classmethod
    # def mod_polarization_log(cls, f, f_max, chirp_mass, z, i, alpha, A_term):
    #     return h

    @classmethod
    def mod_polarization_array(cls, f, f_max, chirp_mass, z, alpha, A_term):
        # returns LOG OF VALUE for plotting purposes
        arr = []

        for i in range(50):
            h_tilde = cmath.log(cls.mod_amplitude(chirp_mass, z, f) *
                                cmath.exp(1j * cls.psi(alpha, A_term, chirp_mass, z, f)[i])) if f < f_max else 0
            arr.append(h_tilde)
        return np.array(arr).real, np.array(arr).imag

    @classmethod
    def chirp_mass_e(cls, chirp_mass, z):
        """ M_e is the source chirp mass, related to measured chirp mass by the following """
        return chirp_mass / (1 + z)

    @classmethod
    def m(cls, chirp_mass, z, m_1, m_2):
        """ calculates symmetric mass ratio eta from source chirp mass and component masses, returns m """
        M_e = cls.chirp_mass_e(chirp_mass, z)
        eta = lambda ma, ss: (ma * ss) / (ma + ss)
        return M_e / (eta(m_1, m_2)**(3/5)), eta(m_1, m_2)

    @classmethod
    def v_lso(cls):
        """ proportional to velocity at last stable orbit """
        return 0

    @classmethod
    def psi_gr_numerical_coefficients(cls, chirp_mass, f, z, m_1, m_2):
        v = cls.u(chirp_mass, f)
        v_lso - cls.v_lso()
        m, eta = m(chirp_mass, z, m_1, m_2)
        gamma = .577216  # euler mascheroni constant (dimensionless)

        v_2 = 20 / 9 * (743 / 336 + (11 / 4 * eta))
        v_3 = 16 * math.pi
        v_4 = 10 * (3058673 / 1016064 + (5329 / 1008 * eta) + (617 / 144 * eta**2))
        v_5 = math.pi * (38645 / 756 - (65 / 9 * eta)) * (1 + 3 * math.ln(v / v_lso))
        v_6 = (11583231236531 / 4694215680
               - (640 / 3 * math.pi**2)
               - (6848 * gamma / 21)
               - (6848 / 21 * math.ln(4 * v))
               + (-15737765635 / 3048192 + 2255 * (math.pi ** 2) / 12) * eta
               + (76055/1728 * eta**2)
               - (127825/1296 * eta**3))
        v_7 = math.pi * (77096675 / 254016 + 378515 / 1512 * eta - 74045 / 756 * eta**2)
        return (v_2 * v**2) - (v_3 * v**3) + (v_4 * v**4) + (v_5 * v**5) + (v_6 * v**6) + (v_7 * v**7)

    """ do not calculate f_dot and set equal to zero to recover max value attained by f
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
