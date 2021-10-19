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
from pycbc.types import FrequencySeries
import pycbc.waveform


class TimeDomain:

    def __init__(self):
        pass

    @classmethod
    def test_waveform(**args
        mass_1 = args['mass1']
        mass_2 = args['mass2']
        df = args['delta_f']
        flow = args['f_lower']
        f = np.linspace(10e-5, 10e-1)
        # t = f / f.max() * (tpeak - tlow) + tlow
        offset = - len(f) * df
        MP = ModifiedPolarization()
        chirp_mass = MP.chirp_mass(mass_1, mass_2)

        wf = MP.std_polarization_array(f, 4, np.linspace(0, 4), chirp_mass, mass_1, mass_2)
        wf_real = FrequencySeries(wf[0], delta_f=df, epoch=offset)
        wf_imag = FrequencySeries(wf[1], delta_f=df, epoch=offset)
        return wf_real, wf_imag

    @classmethod
    def register_standard_waveform(cls):
        approximant = 'standard'
        pycbc.waveform.add_custom_waveform('test', cls.test_waveform, 'frequency', force=True)
        hp, hc = pycbc.waveform.get_fd_waveform(approximant='test', mass1=65, mass2=80, delta_f=1.0 / 4, f_lower=40)

        plt.plot(f, hc)  # plot imag values
        plt.xlabel('$lg(f)$')
        plt.ylabel('$lg(h)$')
        plt.title('Standard polarization in frequency space')
        plt.xlabel('Frequency (Hz)')
        plt.show()
        return cls.plot_pycbc_ifft(approximant, hc)

    @classmethod
    def plot_TaylorF2(cls):
        approximant = 'TaylorF2'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=6, mass2=6,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)
        return cls.plot_ifft_of_waveform(approximant, hp)

    @classmethod
    def plot_IMRPhenomA(cls):
        approximant = 'IMRPhenomA'
        hp, _ = waveform.get_fd_waveform(approximant=approximant,
                                         mass1=65, mass2=80,  # m1 / m2 <= 4 and 50 <= M_e / M_sol <= 200
                                         delta_f=1.0 / 4, f_lower=40)
        return cls.plot_ifft_of_waveform(approximant, hp)

    @classmethod
    def plot_pycbc_ifft(cls, approximant, hp):
        """ gets the frequency domain waveform for a signal and perform ifft to plot over time domain """
        delta_t = 1.0 / 4096
        delta_f = 1.0 / 40

        # perform ifft
        t_len = int(1.0 / delta_f / delta_t)
        # hp = np.resize(hp, (int(t_len / 2 + 1))) # CANT USE NP RESIZE BC THEN ITS NOT A PYCBC ARRAY
        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(hp, sp)

        # plot strain in time domain
        plt.plot(sp.sample_times, sp, label=approximant + ' IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt.show()
