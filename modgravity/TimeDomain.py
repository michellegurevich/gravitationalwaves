import os
import platform
import sys
import camb
import numpy as np
import scipy.fft
from matplotlib import pyplot as plt
from scipy.signal import hilbert
from pycbc import types, fft, waveform
from pycbc.types import FrequencySeries
from ModifiedPolarization import ModifiedPolarization


class TimeDomain:
    # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
    m_1 = 30 * 4.925 * 10e-6
    m_2 = 30 * 4.925 * 10e-6
    df = 1 / 320
    f = np.linspace(30, 350, int(320 / df))
    z_max = 1

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
    def register_test_waveform(cls):
        """ register test waveform to pycbc dictionary """
        return waveform.add_custom_waveform('test', cls.test_waveform, 'frequency', force=True)

    @classmethod
    def get_fd(cls, d):
        hp, hc = waveform.get_fd_waveform(approximant=d['approximant'],
                                         mass1=d['mass1'], mass2=d['mass2'],
                                         delta_f=d['delta_f'], f_lower=d['f_lower'])
        return hp, hc

    @classmethod
    def plot_fd(cls, hp, hc, d):
        plt.plot(hp.sample_frequencies, hp, label=d['approximant'])
        # plt.plot(hc.sample_frequencies, hc, label=d['approximant'])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Strain')
        plt.legend()
        return plt

    @classmethod
    def get_td(cls, d):
        sp, sc = waveform.get_td_waveform(approximant=d['approximant'],
                                          mass1=d['mass1'], mass2=d['mass2'],
                                          delta_t=d['delta_t'], f_lower=d['f_lower'])
        return sp, sc

    @classmethod
    def plot_td(cls, sp, sc, d):
        plt.plot(sp.sample_times, np.abs(sp), label=d['approximant'])
        # plt.plot(sc.sample_times, np.abs(sc), label=d['approximant'])
        plt.xlabel('Time (s)')
        plt.ylabel('Strain')
        plt.legend()
        return plt

    @classmethod
    def plot_pycbc_ifft(cls, approximant, sptilde):
        """ gets the frequency domain waveform for a signal and computes ifft to plot over time domain """
        delta_t = 1.0 / 4096
        delta_f = cls.df

        # perform ifft
        t_len = int(1.0 / delta_t / delta_f)
        # sptilde.resize(t_len / 2 + 1)               # COMMENTING THIS OUT BREAKS TAYLOR F2 BUT NOT TEST

        # generate empty array of size t_len
        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        # use pycbc ifft to get sample times
        fft.ifft(sptilde, sp)
        # use scipy ifft to compute results
        ifft = scipy.fft.irfft(sptilde.data)

        # plot strain in time domain
        plt.plot(sp.sample_times, np.abs(ifft), label=approximant+' (IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt

    """
    
    define waveform dicts and pass these through as args to a single fctn which has a switch statement to match 
    relevant keywords with parameter values
    
    @classmethod
    def plot_TaylorF2(cls):
        approximant = 'TaylorF2'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=cls.df, f_lower=40)
        return cls.plot_pycbc_ifft(approximant, hp)

    @classmethod
    def plot_IMRPhenomA(cls):
        approximant = 'IMRPhenomA'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=65, mass2=80,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)
        return cls.plot_pycbc_ifft(approximant, hp)
        
    """
