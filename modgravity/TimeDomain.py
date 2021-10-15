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


class TimeDomain:

    def __init__(self):
        pass

    @staticmethod
    def plot_modified_waveform_ifft():
        return 0

    @staticmethod
    def test_waveform():
        flow = args['f_lower']  # Required parameter
        df = args['delta_f']  # Required parameter
        fpeak = args['fpeak']  # A new parameter for my model

        t = numpy.arange(0, 10, df)
        f = t / t.max() * (fpeak - flow) + flow
        a = t

        wf = numpy.exp(2.0j * numpy.pi * f * t) * a

        # Return product should be a pycbc time series in this case for
        # each GW polarization
        #
        #
        # Note that by convention, the time at 0 is a fiducial reference.
        # For CBC waveforms, this would be set to where the merger occurs
        offset = - len(t) * dt
        wf = FrequencySeriesSeries(wf, delta_f=df, epoch=offset)
        return wf.real(), wf.imag()

    @staticmethod
    def add_test():
        # This tells pycbc about our new waveform so we can call it from standard
        # pycbc functions. If this were a frequency-domain model, select 'frequency'
        # instead of 'time' to this function call.
        pycbc.waveform.add_custom_waveform('test', test_waveform, 'time', force=True)

        # Let's plot what our new waveform looks like
        hp, hc = pycbc.waveform.get_td_waveform(approximant="test",
                                                f_lower=20, fpeak=50,
                                                delta_t=1.0 / 4096)
        pylab.figure(0)
        pylab.plot(hp.sample_times, hp)
        pylab.xlabel('Time (s)')

        pylab.figure(1)
        hf = hp.to_frequencyseries()
        pylab.plot(hf.sample_frequencies, hf.real())
        pylab.xlabel('Frequency (Hz)')
        pylab.xscale('log')
        pylab.xlim(20, 100)
        return pylab.show()

    @classmethod
    def plot_standard_waveform_ifft(cls):
        MP = ModifiedPolarization()
        f = np.linspace(10e-5, 10e-1)
        z = np.linspace(0, 4)
        z_max = 4
        m_1 = 6
        m_2 = 6
        chirp_mass = MP.chirp_mass(m_1, m_2)
        approximant = 'Waveform'
        hp = MP.std_polarization_array(f, z_max, z, chirp_mass, m_1, m_2)
        return cls.plot_ifft_of_waveform(approximant, hp)

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

        # perform ifft
        t_len = int(1.0 / delta_t / hp.delta_f)
        hp.resize(t_len / 2 + 1)
        sp = types.TimeSeries(types.zeros(t_len), delta_t=delta_t)
        fft.ifft(hp, sp)

        # plot strain in time domain
        plt.plot(sp.sample_times, sp, label=approximant + ' IFFT)')
        plt.ylabel('Strain')
        plt.xlabel('Time (s)')
        plt.legend()
        return plt.show()
