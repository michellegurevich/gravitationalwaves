import camb
import matplotlib
import numpy as np
import os
import platform
import sys
from ModifiedPolarization import ModifiedPolarization
from camb import model, initialpower
from matplotlib import pyplot as plt
from pycbc import types, fft, waveform
from pycbc.waveform import get_td_waveform, td_approximants, fd_approximants
from scipy.signal import hilbert


class TimeDomain:

    def __init__(self):
        pass

    @staticmethod
    def plot_modified_waveform_ifft():
        return 0

    @staticmethod
    def plot_standard_waveform_ifft():
        delta_t = 1.0 / 4096
        z = np.linspace(0, 4)
        m1 = 6
        m2 = 6

        MP = ModifiedPolarization()
        approximant = MP.std_polarization_array(10e-3, 10e-2, 12, z, m1, m2)

        # get TaylorF2 in frequency domain
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=m1, mass2=m2,  # m1 / m2 <= 4 but m1, m2 may be too small for 3.5PN approx
                                         delta_f=1.0 / 4, f_lower=40)

        # perform ifft
        t_len = int(1.0 / delta_t / delta_f)
        hp.resize(t_len / 2 + 1)
        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(hp, sp)

        # plot strain
        plt.plot(sp.sample_times, sp, label='Standard polarization (IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt.show()

    @staticmethod
    def plot_TaylorF2():
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

    @staticmethod
    def plot_IMRPhenomA():
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
