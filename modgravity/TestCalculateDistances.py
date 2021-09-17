import unittest
import numpy as np
from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances


class TestCalculateDistances(unittest.TestCase):
    SC = SetCosmology()
    CD = CalculateDistances()
    z_max = 4
    z = np.linspace(0, z_max)
    z_prime = np.linspace(0, z_max).reshape(-1, 1)
    h0 = 72  # Hubble parameter today in km/s/Mpc
    omega_m = .3
    omega_lambda = .7
    t_e = 5
    t_a = 10

    def test_lum_dist(self):
        self.assertEqual(self.SC.get_ld(self.z), self.CD.lum_dist(self.z_max, self.z_prime, self.h0, self.omega_m,
                                                                  self.omega_lambda))

    def test_lum_dist_integrand(self):
        self.assertEqual(self.SC.get_ld(self.z), self.CD.lum_dist(self.z_max, self.z_prime, self.h0, self.omega_m,
                                                                  self.omega_lambda))

    def test_chi(self):
        # m_g = .0001, E_e = .001, and modification terms are both 0
        self.assertEqual(self.CD.chi(self.z_max, .0001, .001, self.t_e, self.t_a, 0, 0), 154.97909219909434)

if __name__ == '__main__':
    unittest.main()
