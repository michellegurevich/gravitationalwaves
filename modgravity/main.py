import numpy as np
import math
from matplotlib import pyplot as plt
from figures import figures


def main():

    cosmo_params = {
        'redshift': np.linspace(0, 4),
        'frequency': np.linspace(10e0, 10e4),
        'f_e': 1,
        't_e': 1,
        't_a': 2,
    }
    phenom_params = {
        'A': 4,
        'alpha': 4,
        'lambda_g': 1.6 * 10 ** 16,
    }
    wf_params = {
        'approximant': 'TaylorF2',
        'mass1': 10,
        'mass2': 15,
        'delta_f': 1.0 / 170,
        'f_lower': 40,
    }

    f = figures()
    f.phase_check()


if __name__ == '__main__':
    main()
