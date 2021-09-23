import math
import cmath
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
        return math.sqrt(math.pi/30) * (chirp_mass ** 2 / D_L)

    @classmethod
    def mod_amplitude(cls, chirp_mass, z, f):
        """ tilde A(f) term """
        return cls.epsilon * cls.curved_A(chirp_mass, z) * (cls.u(chirp_mass, f) ** (-7/6))

    @classmethod
    def psi_gr(cls, f, chirp_mass):  # , t_c, phi_c):
        freq_term = 2 * math.pi * f * cls.t_c
        numerical_term = 3 / 128 * (cls.u(chirp_mass, f) ** -5/3)  # ADD CROSS PRODUCT
        return freq_term - cls.phi_c - math.pi/4 + numerical_term

    @classmethod
    def psi(cls, alpha, A_term, chirp_mass, z, f):
        psi_gr = cls.psi_gr(f, chirp_mass)  # , t_c, phi_c)
        delta_psi = cls.delta_psi(alpha, A_term, chirp_mass, z, f)
        return psi_gr + delta_psi

    @classmethod
    def mod_polarization(cls, f, f_max, chirp_mass, z, alpha, A_term):
        # needs test to check f > 0
        return cls.mod_amplitude(chirp_mass, z, f) * \
               cmath.exp(1j * cls.psi(alpha, A_term, chirp_mass, z, f)[len(z)-1]) if f < f_max else 0
