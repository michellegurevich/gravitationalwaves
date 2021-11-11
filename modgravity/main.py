import numpy as np
import math
from matplotlib import pyplot as plt

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from ModifiedPolarization import ModifiedPolarization
from TimeDomain import TimeDomain
from Plots import Plots


def main():
    TD = TimeDomain()

    TaylorF2 = {
        'approximant': 'TaylorF2',
        'mass1': 20, 'mass2': 20,
        'delta_f': 1 / 320, 'f_lower': 40,
    }

    TaylorT2 = {
        'approximant': 'TaylorT2',
        'mass1': 20, 'mass2': 20,
        'delta_t': 1.0 / 4096,
        'f_lower': 40
    }

    IMRPhenomA = {
        'approximant': 'IMRPhenomA',
        'mass1': 65, 'mass2': 80,
        'delta_f': 1 / 4, 'f_lower': 40,
        'delta_t': 1.0 / 4096
    }

    test = {
        'approximant': 'test',
        'mass1': TD.m_1, 'mass2': TD.m_2,
        'delta_f': 1 / 320, 'f_lower': 40,
        'delta_t': 1.0 / 4096
    }

    # plot frequency of TaylorF2
    plt.subplot(2, 2, 1)
    hp, hc = TD.get_fd(TaylorF2)
    TD.plot_fd(hp, hc, TaylorF2)

    # plot ifft of TaylorF2
    plt.subplot(2, 2, 2)
    sp, sc = TD.get_td(TaylorT2)
    TD.plot_td(sp, sc, TaylorT2)

    # plot strain against frequency
    plt.subplot(2, 2, 3)
    TD.register_test_waveform()
    kp, kc = TD.get_fd(test)
    TD.plot_fd(kp, kc, test)

    # perform ifft to plot strain against time
    plt.subplot(2, 2, 4)
    TD.plot_pycbc_ifft('test', kp)

    plt.show()


if __name__ == '__main__':
    main()
