import os
import platform
import sys
import camb
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import hilbert
from pycbc import types, fft, waveform
from pycbc.types import FrequencySeries
from ModifiedPolarization import ModifiedPolarization


class TimeDomain:


    def __init__(self):
        pass

    @classmethod
    def test_waveform(cls, **args):
        mass_1 = args['mass1']
        mass_2 = args['mass2']
        df = args['delta_f']
        flow = args['f_lower']
        f = np.linspace(30, 100, int(70 / df))
        # t = f / f.max() * (tpeak - tlow) + tlow
        offset = - len(f) * df
        MP = ModifiedPolarization()

        wf = MP.std_polarization_array(f, 4, np.linspace(0, 4), mass_1, mass_2)
        wf_real = FrequencySeries(wf[0], delta_f=df, epoch=offset)
        wf_imag = FrequencySeries(wf[1], delta_f=df, epoch=offset)
        return wf_real, wf_imag

    @classmethod
    def register_test_waveform(cls):
        approximant = 'test'
        return waveform.add_custom_waveform(approximant, cls.test_waveform, 'frequency', force=True)

    @classmethod
    def get_test_waveform(cls):
        approximant = 'test'
        hp, hc = waveform.get_fd_waveform(approximant=approximant, mass1=65, mass2=80, delta_f=1.0 / 70, f_lower=40)
        return hp, hc

    @classmethod
    def plot_test_waveform(cls, hp, hc, f):
        plt.plot(f, hc)  # plot imag values
        plt.xlabel('$lg(f)$')
        plt.ylabel('$lg(h)$')
        plt.title('Standard polarization in frequency space')
        plt.xlabel('Frequency (Hz)')
        return plt.show()

    @classmethod
    def plot_TaylorF2(cls):
        approximant = 'TaylorF2'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)
        return cls.plot_pycbc_ifft(approximant, hp)

    @classmethod
    def plot_TaylorF2_freq(cls):
        approximant = 'TaylorF2'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)
        plt.plot(hp.sample_frequencies, hp)
        plt.ylabel('Strain')
        plt.xlabel('Frequency (Hz)')
        return plt.show()

    @classmethod
    def plot_TaylorF2(cls):
        approximant = 'TaylorF2'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)
        return cls.plot_pycbc_ifft(approximant, hp)

    @classmethod
    def plot_IMRPhenomA(cls):
        approximant = 'IMRPhenomA'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=65, mass2=80,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)
        return cls.plot_pycbc_ifft(approximant, hp)

    @classmethod
    def plot_pycbc_ifft(cls, approximant, hp):
        """ gets the frequency domain waveform for a signal and perform ifft to plot over time domain """
        delta_t = 1.0 / 4096
        delta_f = 1.0 / 70

        # perform ifft
        t_len = int(1.0 / delta_t / delta_f)
        hp.resize(t_len/2 + 1)

        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(hp, sp)

        # plot strain in time domain
        plt.plot(sp.sample_times, sp, label=approximant + ' IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt.show()
