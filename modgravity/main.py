import numpy as np
from pycbc.types import FrequencySeries
from matplotlib import pyplot as plt
import pycbc.waveform
from ModifiedPolarization import ModifiedPolarization
from Plots import Plots
from TimeDomain import TimeDomain


def test_waveform(**args):
    mass_1 = args['mass1']
    mass_2 = args['mass2']
    df = args['delta_f']
    flow = args['f_lower']

    f = np.linspace(10e-5, 10e-1)
    t = np.linspace(-10, 10)
    # t = f / f.max() * (tpeak - tlow) + tlow

    MP = ModifiedPolarization()
    chirp_mass = MP.chirp_mass(6, 6)
    wf = MP.std_polarization_array(f, 4, np.linspace(0, 4), chirp_mass, 6, 6)
    # wf = MP.mod_polarization_array(f, 10e-1, 4, np.linspace(0, 4), 0.001, 0.01, chirp_mass, 6, 6)

    offset = - len(f) * df
    wf_real = FrequencySeries(wf[0], delta_f=df, epoch=offset)
    wf_imag = FrequencySeries(wf[1], delta_f=df, epoch=offset)
    return wf_real, wf_imag

def main():
    pycbc.waveform.add_custom_waveform('test', test_waveform, 'frequency', force=True)

    hp, hc = pycbc.waveform.get_fd_waveform(approximant='test', mass1=65, mass2=80, delta_f=1.0 / 4, f_lower=40)
    f = np.linspace(10e-5, 10e-1)
    plt.figure(0)
    #plt.plot(f, hc)  # plot imag values
    #plt.xlabel('$lg(f)$')
    #plt.ylabel('$lg(h)$')
    #plt.title('Standard polarization in frequency space')
    #plt.xlabel('Frequency (Hz)')
    #plt.show()

    #plt.figure(1)
    P = Plots()
    MP = ModifiedPolarization()
    #chirp_mass = MP.chirp_mass(6, 6)
    # P.modified_polarization(f, 10e-1, 4, np.linspace(0, 4), 0.001, 0.01, chirp_mass, 6, 6)
    #P.standard_polarization(f, 4, np.linspace(0, 4), chirp_mass, 6, 6)

    TD = TimeDomain()
    TD.plot_pycbc_ifft('test', hc)



if __name__ == '__main__':
    main()
