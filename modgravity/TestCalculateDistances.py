import unittest
import numpy as np
import math

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances


class TestCalculateDistances(unittest.TestCase):
    SC = SetCosmology()
    CD = CalculateDistances()

    z_max = 4
    z = np.linspace(0, z_max)

    alpha = 0
    A_term = 0
    m_g = .0001
    E_e = .001
    t_e = 5
    t_a = 10

    def test_lum_dist(self):
        """ iterates over luminosity distance np arrays (for a given redshift) from camb and my code; returns true if
        each set of corresponding elements (i.e. each set of luminosity distances) agrees within a given tolerance """
        mine = self.CD.lum_dist_array(self.z_max, self.z)
        cambs = self.SC.get_ld(self.z)
        for element in range(len(mine)):
            self.assertTrue(math.isclose(mine[element], cambs[element], rel_tol=.1))

    def test_chi(self):
        self.assertEqual(self.CD.chi(self.z_max, self.alpha, self.A_term, self.m_g, self.E_e, self.t_e, self.t_a), 154.97909219909434)


if __name__ == '__main__':
    unittest.main()
