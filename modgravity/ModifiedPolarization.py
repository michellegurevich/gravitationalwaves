import math
from CalculateDistances import CalculateDistances

class ModifiedPolarization:

    def __init__(self):
        self.h = 1  # Planck constant
        self.f_e = 10 ** 4  # Hz
        self.lambda_g = 1.6 * 10 ** 16  # m

        self.E_e = h * f_e
        self.m_g = h / lambda_g

        self.h0 = 72
        self.omega_m = .3
        self.omega_lambda = .7

    @classmethod
    def lambda_A_term(cls, A_term, alpha):
        return h * A_term ** (1 / (alpha - 2))c # always has units of length irrespective of alpha value

    @classmethod
    def u(cls, chirp_mass, f):
        return math.pi * chirp_mass * f

    @classmethod
    def beta(cls, chirp_mass, z):
        CD = CalculateDistances()
        D_0 = CD.alpha_dist(cls, z, 0, cls.h0, cls.omega_m, cls.omega_lambda)
        return (math.pi ** 2 * D_0 * chirp_mass) / (cls.lambda_g ** 2 (1 + z))
