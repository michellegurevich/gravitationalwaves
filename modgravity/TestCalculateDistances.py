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

    def test_lum_dist(self):
        self.assertEqual(self.SC.get_ld(self.z), self.CD.lum_dist(self.z_max, self.z_prime, self.h0, self.omega_m,
                                                                  self.omega_lambda))


if __name__ == '__main__':
    unittest.main()