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
    def plot_waveform(cls):
        delta_t = 1.0 / 4096
        delta_f = 1.0 / 4

        # Get a frequency domain waveform
        sptilde, sctilde = waveform.get_fd_waveform(approximant="TaylorF2",
                                                    mass1=6, mass2=6, delta_f=1.0 / 4, f_lower=40)

        MP = ModifiedPolarization()
        h_tilde = MP.mod_polarization_array(f, f_max, chirp_mass, z, alpha, A_term)

        # FFT it to the time-domain
        tlen = int(1.0 / delta_t / delta_f)
        h_tilde.resize(tlen / 2 + 1)
        sp = types.TimeSeries(types.zeros(tlen), delta_t=delta_t)
        fft.ifft(h_tilde, h_tilde)  # frequency domain -> time domain

        plt.plot(sp.sample_times, sp, label="TaylorF2 (IFFT)")

        plt.title('FFT of waveform')
        plt.ylabel('$h_{+,x}(t)$')
        plt.xlabel('Time (s)')

        # type(sptilde): pycbc.types.frequencyseries.FrequencySeries models a frequency series consisting of uniformly sampled scalar values
        return plt.show()

    @classmethod
    def plot_phenomA(cls):
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
