"""Render the SPA-specific examples relocated from former image chapters."""
from pathlib import Path
import sys,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import skimage as ski
import mrcfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'book'))
from _support.spa_filters import ctf_1d,radial_frequency_grid,phase_flip_spectrum,wiener_ctf_spectrum
out=ROOT/'book/images/spa';s=np.linspace(0,.5,4000)
fig,ax=plt.subplots(figsize=(9,3.5))
ax.plot(s,ctf_1d(s,15000),label='Pure phase, no envelope')
ax.plot(s,ctf_1d(s,15000,amplitude_contrast=.1,phase_shift_rad=.15,b_factor_A2=40),label='Amplitude contrast + phase shift + envelope')
ax.axhline(0,color='gray',lw=.6);ax.set(xlabel='Spatial frequency (Å⁻¹)',ylabel='CTF');ax.legend();fig.tight_layout();fig.savefig(out/'ctf_terms_migrated.png');plt.close(fig)
with mrcfile.mmap(ROOT/'book/data/70S_Conform1.mrc',permissive=True,mode='r') as m:
    projection=np.asarray(m.data,dtype=float).sum(axis=0);pixel=float(m.voxel_size.x);shape=list(m.data.shape)
projection=(projection-projection.min())/(projection.max()-projection.min())
H=ctf_1d(radial_frequency_grid(projection.shape,pixel),15000,amplitude_contrast=.1,b_factor_A2=30)
rng=np.random.default_rng(7);noise=rng.normal(scale=.03*projection.std(),size=projection.shape)
observed=np.fft.ifft2(H*np.fft.fft2(projection)).real+noise
F=np.fft.fft2(observed);flip=np.fft.ifft2(phase_flip_spectrum(F,H)).real;wiener=np.fft.ifft2(wiener_ctf_spectrum(F,H,.03)).real
fig,axes=plt.subplots(1,4,figsize=(12,3))
for ax,a,title in zip(axes,[projection,observed,flip,wiener],['70S projection','CTF + noise','Phase flipping','Wiener-style']):
    ax.imshow(a,cmap='gray');ax.set_title(title);ax.axis('off')
fig.tight_layout();fig.savefig(out/'ctf_correction_migrated.png');plt.close(fig)
phantom=ski.transform.resize(ski.data.shepp_logan_phantom(),(96,128),anti_aliasing=True)
f1=np.fft.fft(phantom.sum(axis=0));f2=np.fft.fft2(phantom)[0,:]
assert np.allclose(f1,f2,atol=1e-10)
fig,axes=plt.subplots(1,2,figsize=(9,3.5));axes[0].imshow(phantom,cmap='gray');axes[0].set_title('2D object');axes[0].axis('off');axes[1].plot(np.abs(np.fft.fftshift(f1)),label='FFT of projection');axes[1].plot(np.abs(np.fft.fftshift(f2)),'--',label='Central slice');axes[1].legend();fig.tight_layout();fig.savefig(out/'slice_example_migrated.png');plt.close(fig)
report={'source':'migrated from pre-restoration 03_fourier; CTF demonstration now uses real 70S forward projection','seed':7,'volume_shape':shape,'pixel_size_A':pixel,'df_A':15000,'voltage_kv':300,'amplitude_contrast':.1,'B_A2':30,'noise_sd_fraction':.03,'K':.03,'slice_max_error':float(np.max(np.abs(f1-f2))),'phase_flip_magnitude_error':float(np.max(np.abs(np.abs(phase_flip_spectrum(F,H))-np.abs(F))))}
(ROOT/'notes/image_revision_20260913/spa_migration.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
