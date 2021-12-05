import unittest
import numpy as np
import math

from ModifiedPolarization import ModifiedPolarization


class TestModel(unittest.TestCase):
    MP = ModifiedPolarization()

    def test_lambda_A_term(self):
        alpha = 4
        A_term = .001
        actual = MP.lambda_A_term(alpha, A_term)
        expected = 1 * A_term ** (1 / (alpha-2))
        self.assertEqual(actual, expected)

    def test_u(self):
        actual = MP.u(f, chirp_mass)
        expected = 0
        self.assertEqual(actual, expected)

    def test_beta(self):
        actual = MP.beta(z_max, z, chirp_mass)
        expected = 0
        self.assertEqual(actual, expected)

    def test_zeta_if(self):
        """ case where alpha = 1 """
        actual = MP.zeta(z_max, z, 1, A_term, chirp_mass)
        expected = 0
        self.assertEqual(actual, expected)

    def test_zeta_else(self):
        """ case where alpha != 1 """
        actual = MP.zeta(z_max, z, 2, A_term, chirp_mass)
        expected = 0
        self.assertEqual(actual, expected)

    def test_chirp_mass(self):
        actual = MP.chirp_mass(m_1, m_2)
        expected = 0
        self.assertEqual(actual, expected)

    def test_chirp_mass_e(self):
        actual = MP.chirp_mass_e(z_max, chirp_mass)
        expected = 0
        self.assertEqual(actual, expected)

    def test_m_sum(self):
        actual, _ = MP.m(m_1, m_2)
        expected = 0
        self.assertEqual(actual, expected)

    def test_eta(self):
        _, actual = MP.m(m_1, m_2)
        expected = 0
        self.assertEqual(actual, expected)

    def test_v_lso(self):
        actual = MP.v_lso()
        expected = 1 / math.sqrt(6)
        self.assertEqual(actual, expected)

    def test_curved_A(self):
        actual = MP.curved_A(z_max, x, chirp_mass)
        expected = 0
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
