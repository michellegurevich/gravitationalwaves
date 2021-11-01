import numpy as np
import math
from matplotlib import pyplot as plt

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from ModifiedPolarization import ModifiedPolarization
from TimeDomain import TimeDomain
from Plots import Plots


def main():
    z_max = 4
    z = np.linspace(0, z_max)
    # alpha = 3
    # A_term = .0001
    delta_f = 1/70
    f = np.linspace(30, 100, int(70 / delta_f))
    # f_cut = 10e-1
    m_1 = 30 * 4.925 * 10e-6
    m_2 = 30 * 4.925 * 10e-6

    P = Plots()
    # P.modified_polarization(f, f_cut, z_max, z, alpha, A_term, chirp_mass, m_1, m_2)
    P.standard_polarization(f, z_max, z, m_1, m_2)

    TD = TimeDomain()

    # plot frequency of TaylorF2
    TD.plot_TaylorF2_freq()

    # plot ifft of TaylorF2
    TD.plot_TaylorF2()

    # register test waveform to pycbc fd collection
    TD.register_test_waveform()

    # get waveform from pycbc fd collection
    hp, hc = TD.get_test_waveform()

    # plot strain against frequency
    TD.plot_test_waveform(hp, hc, f)

    # perform ifft to plot strain against time
    TD.plot_pycbc_ifft('test', hp)


if __name__ == '__main__':
    main()
