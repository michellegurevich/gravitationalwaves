import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower

from pycbc.waveform import get_td_waveform, td_approximants, fd_approximants
from pycbc import types, fft, waveform
from scipy.signal import hilbert

from ModifiedPolarization import ModifiedPolarization


class TimeDomain:

    def __init__(self):
        pass

    @classmethod
    def plot_TaylorF2(cls):
        delta_t = 1.0 / 4096
        delta_f = 1.0 / 4

        # get TaylorF2 in frequency domain
        hp, _ = waveform.get_fd_waveform(approximant="TaylorF2",
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 but m1, m2 may be too small for 3.5PN approx
                                         delta_f=1.0 / 4, f_lower=40)

        # perform ifft
        t_len = int(1.0 / delta_t / delta_f)
        hp.resize(t_len / 2 + 1)
        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(hp, sp)

        # plot strain
        plt.plot(sp.sample_times, sp, label='TaylorF2 (IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt.show()

    @classmethod
    def plot_IMRPhenomA(cls):
        """ gets the frequency domain waveform for phenomA signa and iffts it to plot over time domain """
        delta_t = 1.0 / 4096
        delta_f = 1.0 / 4

        # get phenomA in frequency domain
        hp, _ = waveform.get_fd_waveform(approximant='IMRPhenomA',
                                         mass1=65, mass2=80,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)

        # perform ifft
        t_len = int(1.0 / delta_t / delta_f)
        hp.resize(t_len / 2 + 1)
        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(hp, sp)

        # plot strain
        plt.plot(sp.sample_times, sp, label='IMRPhenomA (IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt.show()
