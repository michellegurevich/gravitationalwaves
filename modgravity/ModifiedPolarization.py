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
    t_c = .0001  # verify in references
    phi_c = .0001  # verify in references

    h0 = 72
    omega_m = .3
    omega_lambda = .7

    def __init__(self):
        pass

    @classmethod
    def lambda_A_term(cls, alpha, A_term):
        return cls.h * A_term ** (1 / (alpha - 2))  # always has units of length irrespective of alpha value

    @classmethod
    def u(cls, f, chirp_mass):
        return math.pi * chirp_mass * f

    @classmethod
    def beta(cls, z_max, z, chirp_mass):
        CD = CalculateDistances()
        D_0 = CD.alpha_dist(z_max, z, 0)
        return (math.pi ** 2 * D_0 * chirp_mass) / (cls.lambda_g ** 2 * (1 + z_max))

    @classmethod
    def zeta(cls, z_max, z, alpha, A_term, chirp_mass):
        CD = CalculateDistances()
        lambda_A_term = cls.lambda_A_term(alpha, A_term)
        if alpha == 1:
            D_1 = CD.alpha_dist(z_max, z, alpha)
            zeta = (math.pi * D_1) / lambda_A_term
        else:
            D_alpha = CD.alpha_dist(z_max, z, alpha)
            term_i = math.pi ** (2 - alpha) / (1 - alpha)
            term_ii = D_alpha / (lambda_A_term ** (2 - alpha))
            term_iii = (chirp_mass ** (1 - alpha)) / ((1 + z_max) ** (1 - alpha))
            zeta = term_i * term_ii * term_iii
        return zeta

    @classmethod
    def curved_A(cls, z_max, z, chirp_mass):
        CD = CalculateDistances()
        D_L = CD.lum_dist(z_max, z)
        return math.sqrt(math.pi / 30) * (chirp_mass ** 2 / D_L)

    @classmethod
    def mod_amplitude(cls, f, z_max, z, chirp_mass):
        """ A~(f) term """
        return cls.epsilon * cls.curved_A(z_max, z, chirp_mass) * (cls.u(f, chirp_mass) ** (-7 / 6))

    @classmethod
    def delta_psi(cls, f, z_max, z, alpha, A_term, chirp_mass):
        term_i = cls.beta(z_max, z, chirp_mass) / cls.u(f, chirp_mass)
        if alpha == 1:
            term_ii = cls.zeta(z_max, z, alpha, A_term, chirp_mass) * math.log(cls.u(f, chirp_mass))
        else:
            term_ii = cls.zeta(z_max, z, alpha, A_term, chirp_mass) * (cls.u(f, chirp_mass) ** (alpha - 1))
        return -1 * (term_i + term_ii)

    @classmethod
    def psi_gr(cls, f, z_max, chirp_mass, m_1, m_2):
        freq_term = 2 * math.pi * f * cls.t_c
        table_coefficient_terms = cls.psi_gr_numerical_coefficients(f, z_max, chirp_mass, m_1, m_2)
        numerical_term = 3 / 128 * (cls.u(f, chirp_mass) ** -5 / 3) + table_coefficient_terms
        return freq_term - cls.phi_c - math.pi / 4 + numerical_term

    @classmethod
    def psi(cls, f, z_max, z, alpha, A_term, chirp_mass, m_1, m_2):
        psi_gr = cls.psi_gr(f, z_max, chirp_mass, m_1, m_2)
        delta_psi = cls.delta_psi(f, z_max, z, alpha, A_term, chirp_mass)
        return psi_gr + delta_psi

    @classmethod
    def mod_polarization_array(cls, f, f_cut, z_max, z, alpha, A_term, chirp_mass, m_1, m_2):
        """ assigns the LOG OF VALUE (for plotting purposes) of the product of modified amplitude and exp(i*psi) to an
        array whose length corresponds to that of psi """
        arr = []
        for i in range(50):
            h_tilde = np.log(cls.mod_amplitude(f[i], z_max, z, chirp_mass) *
                             np.exp(1j * cls.psi(f[i], z_max, z, alpha, A_term, chirp_mass, m_1, m_2))) \
                if f[i] < f_cut else 0
            arr.append(h_tilde)
        return [arr[i].real for i in range(50)], [arr[i].imag for i in range(50)]

    @classmethod
    def chirp_mass(cls, m_1, m_2):
        """ M can be calculated from component masses m_1 and m_2 by the following """
        return ((m_1 * m_2) ** (3/5)) / ((m_1 + m_2) ** (1/5))

    @classmethod
    def chirp_mass_e(cls, z_max, chirp_mass):
        """ M_e is the (emitted) source chirp mass, related to (observed) measured chirp mass by the following """
        return chirp_mass / (1 + z_max)

    @classmethod
    def m(cls, z_max, chirp_mass, m_1, m_2):
        """ calculates symmetric mass ratio eta from source chirp mass and component masses; returns m and eta """
        M_e = cls.chirp_mass_e(z_max, chirp_mass)
        eta = (m_1 * m_2) / ((m_1 + m_2) ** 2)
        return M_e / (eta**(3/5)), eta

    @classmethod
    def v_lso(cls):
        """ proportional to velocity at last stable orbit in the strongly relativistic regime (here we use the
        Schwarzschild metric termination point) """
        return 1 / math.sqrt(6)

    @classmethod
    def psi_gr_numerical_coefficients(cls, f, z_max, chirp_mass, m_1, m_2):
        v = cls.u(f, chirp_mass)
        v_lso = cls.v_lso()
        _, eta = cls.m(z_max, chirp_mass, m_1, m_2)
        gamma = .577216  # euler mascheroni constant (dimensionless)

        v_2 = 20 / 9 * (743 / 336 + (11 / 4 * eta))
        v_3 = 16 * math.pi
        v_4 = 10 * (3058673 / 1016064 + (5329 / 1008 * eta) + (617 / 144 * eta**2))
        v_5 = math.pi * (38645 / 756 - (65 / 9 * eta)) * (1 + 3 * np.log(v / v_lso))
        v_6 = (11583231236531 / 4694215680
               - (640 / 3 * math.pi**2)
               - (6848 * gamma / 21)
               - (6848 / 21 * np.log(4 * v))
               + (-15737765635 / 3048192 + 2255 * (math.pi ** 2) / 12) * eta
               + (76055/1728 * eta**2)
               - (127825/1296 * eta**3))
        v_7 = math.pi * (77096675 / 254016 + 378515 / 1512 * eta - 74045 / 756 * eta**2)

        return (v_2 * v**2) - (v_3 * v**3) + (v_4 * v**4) + (v_5 * v**5) + (v_6 * v**6) + (v_7 * v**7)

    @classmethod
    def std_polarization_array(cls, f, z_max, z, chirp_mass, m_1, m_2):
        """ assigns the LOG OF VALUE (for plotting purposes) of the product of standard amplitude and exp(i*psi_GR) to
        an array whose length corresponds to that of psi """
        arr = []

        for i in range(50):
            h = np.log(cls.curved_A(z_max, z, chirp_mass) * np.exp(1j * cls.psi_gr(f[i], z_max, chirp_mass, m_1, m_2)))
            arr.append(h)
        return [arr[i].real for i in range(50)], [arr[i].imag for i in range(50)]

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
