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
    m_1 = 10 * 4.925 * 10e-6
    m_2 = 15 * 4.925 * 10e-6
    z_max = 1
    MP = ModifiedPolarization()

    def __init__(self):
        pass

    @classmethod
    def test_waveform(cls, **args):
        df = args['delta_f']
        flow = args['f_lower']

        f = np.linspace(30, 200, int(170 / df))
        offset = - len(f) * df
        wf, wp = cls.MP.std_polarization_array(f, cls.z_max, np.linspace(0, cls.z_max), cls.m_1, cls.m_2)
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
    def plot_pycbc_ifft(cls, d):
        """ gets the frequency domain waveform for a signal and computes ifft to plot over time domain """
        sptilde, sctilde = cls.get_fd(d)

        delta_t = 1.0 / 4096
        delta_f = sctilde.delta_f

        # perform ifft
        t_len = int(1.0 / delta_t / delta_f)
        sctilde.resize(t_len / 2 + 1)

        # use pycbc ifft (with numpy as default backend) to get sample times
        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(sctilde, sp)

        # plot strain in time domain
        plt.plot(sp.sample_times, np.abs(sp), label=d['approximant'] + ' (IFFT)')
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
                                        delta_f=1/320, f_lower=40)
       return cls.plot_pycbc_ifft(approximant, hp)

   @classmethod
   def plot_IMRPhenomA(cls):
       approximant = 'IMRPhenomA'
       hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                        mass1=65, mass2=80,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                        delta_f=1.0 / 4, f_lower=40)
       return cls.plot_pycbc_ifft(approximant, hp)

    # plot frequency of TaylorF2
    plt.subplot(2, 2, 1)
    hp, hc = TD.get_fd(TaylorF2)
    TD.plot_fd(hp, hc, TaylorF2)

    # plot ifft of TaylorF2
    plt.subplot(2, 2, 2)g
    # sp, sc = TD.get_td(TaylorT2)
    # TD.plot_td(sp, sc, TaylorT2)
    TD.plot_pycbc_ifft(TaylorF2)

    # plot strain against frequency
    plt.subplot(2, 2, 3)
    TD.register_test_waveform()
    kp, kc = TD.get_fd(test)
    TD.plot_fd(kp, kc, test)

    # perform ifft to plot strain against time
    plt.subplot(2, 2, 4)
    TD.plot_pycbc_ifft(test)

    plt.show()

    # MP.phase_check()

   """
