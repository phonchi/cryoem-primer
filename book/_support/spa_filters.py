"""SPA-only CTF helpers relocated out of the general image-processing chapters."""
import numpy as np

def electron_wavelength_A(voltage_kv):
    if voltage_kv<=0:raise ValueError('voltage must be positive')
    v=voltage_kv*1000
    return 12.2639/np.sqrt(v*(1+0.97845e-6*v))

def ctf_1d(s,defocus_A,voltage_kv=300.,cs_mm=2.7,amplitude_contrast=0.,phase_shift_rad=0.,b_factor_A2=0.):
    if not 0<=amplitude_contrast<1:raise ValueError('invalid amplitude contrast')
    wavelength=electron_wavelength_A(voltage_kv);s=np.asarray(s)
    gamma=np.pi*wavelength*defocus_A*s*s-.5*np.pi*cs_mm*1e7*wavelength**3*s**4+phase_shift_rad
    return -np.exp(-.25*b_factor_A2*s*s)*(np.sqrt(1-amplitude_contrast**2)*np.sin(gamma)+amplitude_contrast*np.cos(gamma))

def radial_frequency_grid(shape,pixel_size_A):
    y,x=np.meshgrid(np.fft.fftfreq(shape[0],d=pixel_size_A),np.fft.fftfreq(shape[1],d=pixel_size_A),indexing='ij')
    return np.hypot(y,x)

def phase_flip_spectrum(observed,ctf):return np.sign(ctf)*observed

def wiener_ctf_spectrum(observed,ctf,regularization):
    if regularization<=0:raise ValueError('regularization must be positive')
    return np.conj(ctf)*observed/(np.abs(ctf)**2+regularization)
