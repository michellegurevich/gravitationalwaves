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

    # plot frequency of TaylorF2
    plt.subplot(2, 2, 1)
    TD.plot_TaylorF2_freq()

    # plot ifft of TaylorF2
    plt.subplot(2, 2, 2)
    TD.plot_TaylorF2()

    # plot strain against frequency
    plt.subplot(2, 2, 3)
    hp, hc = TD.plot_test_waveform()

    # perform ifft to plot strain against time
    plt.subplot(2, 2, 4)
    TD.plot_pycbc_ifft('test', hp)

    plt.show()


if __name__ == '__main__':
    main()
