import math
from CalculateDistances import CalculateDistances

class ModifiedPolarization:

    def __init__(self):
        self.h = 1  # Planck constant
        self.f_e = 10 ** 4  # Hz
        self.lambda_g = 1.6 * 10 ** 16  # m

        self.E_e = self.h * self.f_e
        self.m_g = self.h / self.lambda_g

        self.h0 = 72
        self.omega_m = .3
        self.omega_lambda = .7

    @classmethod
    def lambda_A_term(cls, A_term, alpha):
        return h * A_term ** (1 / (alpha - 2))  # always has units of length irrespective of alpha value

    @classmethod
    def u(cls, chirp_mass, f):
        return math.pi * chirp_mass * f

    @classmethod
    def beta(cls, chirp_mass, z):
        CD = CalculateDistances()
        D_0 = CD.alpha_dist(cls, z, 0, cls.h0, cls.omega_m, cls.omega_lambda)
        return (math.pi ** 2 * D_0 * chirp_mass) / (cls.lambda_g ** 2 * (1 + z))

    @classmethod
    def zeta(cls, alpha, A_term, chirp_mass):
        lambda_A_term = cls.lambda_A_term(A_term, alpha)
        CD = CalculateDistances()
        if alpha == 1:
            D_1 = CD.alpha_dist(cls, z, 1, cls.h0, cls.omega_m, cls.omega_lambda)
            zeta = (math.pi * D_1) / lambda_A_term
        else:
            D_alpha = CD.alpha_dist(cls, z, alpha, cls.h0, cls.omega_m, cls.omega_lambda)
            term_i = math.pi ** (2 - alpha) / (1 - alpha)
            term_ii = D_alpha / (lambda_A_term ** (2 - alpha))
            term_iii = (chirp_mass ** (1 - alpha)) / ((1 + z) ** (1 - alpha))
            zeta = term_i * term_ii * term_iii
        return zeta

    @classmethod
    def delta_psi(cls, alpha, chirp_mass, z, f):
        term_i = cls.beta(chirp_mass, z) / cls.u(chirp_mass, f)
        if alpha == 1:
            term_ii = cls.zeta(alpha, A_term, chirp_mass) * math.ln(cls.u(chirp_mass, f))
        else:
            term_ii = cls.zeta(alpha, A_term, chirp_mass) * (cls.u(chirp_mass, f) ** (alpha - 1))
        return -1 * (term_i + term_ii)