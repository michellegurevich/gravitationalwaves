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
    m_1 = 30 * 4.925 * 10e-6
    m_2 = 30 * 4.925 * 10e-6
    df = 1 / 320
    f = np.linspace(30, 350, int(320 / df))
    z_max = 4

    def __init__(self):
        pass

    @classmethod
    def test_waveform(cls, **args):
        df = args['delta_f']
        flow = args['f_lower']

        MP = ModifiedPolarization()
        offset = - len(cls.f) * df
        wf, wp = MP.std_polarization_array(cls.f, cls.z_max, np.linspace(0, cls.z_max), cls.m_1, cls.m_2)
        wf = FrequencySeries(wf, delta_f=df, epoch=offset)
        wp = FrequencySeries(wp, delta_f=df, epoch=offset)
        return wf, wp

    @classmethod
    def plot_test_waveform(cls):
        # register test waveform to pycbc dictionary
        waveform.add_custom_waveform('test', cls.test_waveform, 'frequency', force=True)

        # get test waveform from pycbc dictionary
        hp, hc = waveform.get_fd_waveform(approximant='test', delta_f=cls.df, f_lower=40)

        # plot (either real or imag values of) test waveform in frequency space
        #plt.plot(cls.f, hp)  # plot real values
        #plt.xlabel('Frequency (Hz)')
        #plt.ylabel('Strain')
        #return hp, hc
        return cls.plot_pycbc_ifft('test', hp)

    @classmethod
    def plot_TaylorF2(cls):
        approximant = 'TaylorF2'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=cls.df, f_lower=40)
        return cls.plot_pycbc_ifft(approximant, hp)

    @classmethod
    def plot_TaylorF2_freq(cls):
        approximant = 'TaylorF2'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=cls.df, f_lower=40)
        plt.plot(hp.sample_frequencies, hp)
        plt.ylabel('Strain')
        plt.xlabel('Frequency (Hz)')
        return plt

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
        delta_f = cls.df

        # perform ifft
        t_len = int(1.0 / delta_t / delta_f)
        hp.resize(t_len/2 + 1)  # COMMENTING THIS OUT BREAKS TAYLOR F2 BUT NOT TEST

        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(hp, sp)

        # plot strain in time domain
        plt.plot(sp.sample_times, sp, label=approximant+' (IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt
