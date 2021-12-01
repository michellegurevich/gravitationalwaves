import unittest
import numpy as np
import math

from Model import Model


class TestModel(unittest.TestCase):
    M = Model()

    def test_get_waveform(self):
        pycbc_wf = {
            'approximant': 'TaylorF2',
            'mass1': 10, 'mass2': 15,
            'delta_f': 1.0 / 170, 'f_lower': 40,
        }

        self.assertTrue(M.get_waveform() == pycbc_wf)

    #def test_decompose_waveform(self):
    #    self.assertEqual


if __name__ == '__main__':
    unittest.main()
