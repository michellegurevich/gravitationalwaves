#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower


# In[2]:


<<<<<<< HEAD
from pycbc.waveform import get_td_waveform, td_approximants, fd_approximants
from pycbc import waveform
from scipy.signal import hilbert
=======
from pycbc.waveform import get_td_waveform
from pycbc import waveform
>>>>>>> 90f68f0c6bfa5bfcfac7e61ac11397f3854f3c41


# In[3]:


camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))
sys.path.insert(0,camb_path)
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))


# In[4]:


params = camb.CAMBparams()
# set Hubble paramter to present day (z=0) value
params.set_cosmology(H0 = 67.4)
results = camb.get_background(params, no_thermo = False)

print('derived parameter dictionary: %s'%results.get_derived_params()) # check age should be 13.8 billion years today


# In[5]:


z = 4
add = results.angular_diameter_distance(z) # ang diam distance to object at redshift z = 4
ld = results.luminosity_distance(z) # lum distance to same object (uses relationship to add, not apparent and absolute magnitudes)
print("angular diameter distance to object:", add, "luminosity distance to object:", ld)


# In[6]:


# constants
c = 2.99792458e8 # m/s
G = 6.67259e-11 # m^3/kg/s^2 
M_sun = 1.989e30 # kg

## from www.gw-openscience.org
## source parameters
#m1 = 30 * M_sun
#m2 = 25 * M_sun
#M_total = m1 + m2
#eta = (m1*m2) / M_total**2 # symmetric mass ratio
#mchirp = M_total * eta**(3./5.) # chirp mass

## final BH mass is typically 95% of the total initial mass
#M_total_final = 0.95 * M_total
## Final BH radius
#R_final = 2 * G * M_total_final * M_sun / c**2/1000. # km


<<<<<<< HEAD
# In[7]:


# List of td approximants that are available
#print(td_approximants())

# List of fd citapproximants that are currently available
#print(fd_approximants())


# In[8]:
=======
# In[ ]:





# In[7]:
>>>>>>> 90f68f0c6bfa5bfcfac7e61ac11397f3854f3c41


# plot time domain waveform

hp, hc = get_td_waveform(approximant = 'IMRPhenomC', mass1 = 30, mass2 = 25,
                         spin1z = 0.9, # z component of the first binary component’s dimensionless spin
                         delta_t = 1.0 / 4096, # time step used (from 4096 Hz frequency)
                         f_lower = 40 # starting frequency of waveform in Hz
<<<<<<< HEAD
                         #distance = 10 # luminosity distance in Mpc
=======
>>>>>>> 90f68f0c6bfa5bfcfac7e61ac11397f3854f3c41
                        )
plt.plot(hp.sample_times, hp)
    
plt.title('Time domain of coalescence')
plt.ylabel('$h_{+,x}(t)$')
plt.xlabel('Time (s)')
plt.annotate('Inspiral', (-.4,-2e-19))
plt.annotate('Merger', (-.18,-5e-19))
plt.annotate('Ringdown', (.05,-1e-19))
plt.show()


<<<<<<< HEAD
# In[9]:
=======
# In[8]:
>>>>>>> 90f68f0c6bfa5bfcfac7e61ac11397f3854f3c41


# plot frequency evolution of time domain waveform

<<<<<<< HEAD
hp, hc = waveform.get_td_waveform(approximant='SEOBNRv4HM', mass1 = 30, mass2 = 25, phase_order = 7, 
=======
hp, hc = waveform.get_td_waveform(approximant='SEOBNRv2', mass1 = 30, mass2 = 25, phase_order = 7, 
>>>>>>> 90f68f0c6bfa5bfcfac7e61ac11397f3854f3c41
                                  delta_t = 1.0 / 4096, f_lower = 40)

hp, hc = hp.trim_zeros(), hc.trim_zeros()
amp = waveform.utils.amplitude_from_polarizations(hp, hc)
f = waveform.utils.frequency_from_polarizations(hp, hc)

plt.plot(f.sample_times, f)

<<<<<<< HEAD
plt.title('Frequency evolution')
=======
>>>>>>> 90f68f0c6bfa5bfcfac7e61ac11397f3854f3c41
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.show()


<<<<<<< HEAD
# In[10]:


# from scipy docu

t = np.asarray(hp.sample_times)
t0 = t
dt = 1.0 / 4096

analytic_signal = hilbert(np.asarray(hp))
amplitude_envelope = np.abs(analytic_signal)
instantaneous_phase = np.unwrap(np.angle(analytic_signal))
instantaneous_frequency = (np.diff(instantaneous_phase) / (2.0 * np.pi * dt))


# In[11]:


order = 1
fit = np.polyfit(t0, instantaneous_phase, order)
#print(fit)


# In[12]:


p = np.poly1d(fit) # fit constant offset to the instaneous phase
phi_c = fit[0]
# if order > 0: print(fit[1] / 2 / np.pi)
# print(np.mod(phi_c,2*np.pi))
estimated = p(t0) # re-evaluate the offset term using the fitted values
offsetTerm = estimated
demodulated = instantaneous_phase - offsetTerm
signal = np.cos(demodulated + estimated) * np.abs(analytic_signal)


# In[13]:


plt.figure()
plt.plot(t0,signal, label = 'measured signal') #demodulated signal
plt.plot(t0,np.cos(demodulated)*np.abs(analytic_signal), label = 'demodulated signal') 
plt.legend()
plt.title('Demodulated signal')
plt.ylabel('$h_{+,x}(t)$')
plt.xlabel('Time (s)')
plt.show()
=======
# In[ ]:





# In[ ]:



>>>>>>> 90f68f0c6bfa5bfcfac7e61ac11397f3854f3c41

