import unittest
import numpy as np
import math

from set_cosmology import set_cosmology
from distances import distances


COSMO_PARAMS = {
    'redshift' : np.linspace(0,4),
    'frequency': np.linspace(10e-1,10e-4),
    't_e' : 5,
    't_a' : 10,
}
PHENOM_PARAMSt = {
    'A' : 0,
    'alpha' : 0,
    'lambda_g' : 1.6 * 10 ** 16,
}


class TestCalculateDistances(unittest.TestCase):

    def __init__(self):
        self.cosmology = set_cosmology()
        self.dist = distances(COSMO_PARAMS, PHENOM_PARAMS)

    def test_lum_dist(self):
        """ confirms camb and my code agree on luminosity distances calculation
        within a given tolerance """
        mine = self.dist.luminosity_distance()
        cambs = self.cosmology.get_ld(self.z)
        for element in range(len(mine)):
            self.assertTrue(math.isclose(mine[element], cambs[element], rel_tol=.1))

    def test_chi(self):
        self.assertEqual(self.dist.chi_term(self.z), 154.97909219909434)


if __name__ == '__main__':
    unittest.main()
