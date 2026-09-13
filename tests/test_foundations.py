"""Original Fourier and wavelet experiments; SPA helper checks live separately."""
import numpy as np
import pywt

def test_original_basis_and_real_symmetric_masks(restored_chapter):
    c=restored_chapter('03_fourier.py')
    fn=c['coefficient_basis']
    for kx,ky in [(0,0),(5,2),(-32,17),(64,0)]:
        for kind in ['sine','cosine']:
            result=fn(128,kx,ky,.6,kind)
            # helper returns (frequency coefficients, real spatial image).
            H,h=result
            assert np.max(np.abs(np.fft.ifft2(H).imag))<1e-12
            np.testing.assert_allclose(h,np.fft.ifft2(H).real,atol=1e-12)
    for shape in [(80,80),(79,82)]:
        for cutoff in [1,20,24]:
            M=c['square_lowpass'](shape,cutoff)
            np.testing.assert_array_equal(M,M[np.ix_((-np.arange(shape[0]))%shape[0],(-np.arange(shape[1]))%shape[1])])

def test_all_original_cutoffs_and_feature_methods_present(restored_chapter):
    c=restored_chapter('03_fourier.py')
    assert list(c['cutoffs'])==list(range(1,25))
    assert len(c['lowpass_results'])==len(c['highpass_results'])==24
    for lo,hi in zip(c['lowpass_results'],c['highpass_results']):
        np.testing.assert_allclose(lo+hi,c['gengar'],atol=1e-10)
    assert c['harris_response'].ndim==2
    assert all(x.shape[1]==3 for x in c['blobs_list'])

def test_original_ecg_and_lucario_wavelets(restored_chapter):
    c=restored_chapter('04_wavelet.py')
    np.testing.assert_allclose(c['reconstructed_signals'][:c['signals'].size],c['signals'],atol=1e-9)
    assert c['original'].shape==(110,80,3)
    assert all(np.isfinite(c[k]) for k in ['psnr_noisy','psnr_bayes','psnr_visushrink','psnr_visushrink2','psnr_visushrink4'])
    assert c['psnr_bayes']>c['psnr_noisy']

def test_spa_phase_flip_helper_is_not_in_general_chapter():
    from pathlib import Path
    import sys
    book=Path(__file__).resolve().parents[1]/'book';sys.path.insert(0,str(book))
    from _support.spa_filters import phase_flip_spectrum,wiener_ctf_spectrum
    x=np.ones((8,8),complex)*(2+1j);H=np.linspace(-1,1,64).reshape(8,8)
    np.testing.assert_allclose(np.abs(phase_flip_spectrum(x,H)),np.abs(x))
    assert np.isfinite(wiener_ctf_spectrum(x,H,.03)).all()
    assert 'def phase_flip_spectrum' not in (book/'03_fourier.py').read_text()
