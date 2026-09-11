"""Independent numerical checks for the corrected SPA textbook examples."""
from pathlib import Path
import json
import numpy as np

result={}
# Unit conversions, without assuming all transferred energy stays in thin ice.
result['particles_per_square_micron']=2*8e-17/1e6*6.022e23
result['particles_per_hole']=result['particles_per_square_micron']*np.pi*.6**2
result['holes_for_200000']=200000/result['particles_per_hole']
result['dose_upper_scale_Gy']=2.4*1e17*1e6*1.602176634e-19*1000
n=10000; rho=.01
result['correlated_mean_variance_factor']=(1+(n-1)*rho)/n
result['effective_independent_count']=n/(1+(n-1)*rho)
chi=np.pi/2
result['weak_phase_cross_term_sign']=2*np.real(1j*np.exp(-1j*chi))
assert np.isclose(result['weak_phase_cross_term_sign'],2)
wave=.0197; resolution=3.5
result['first_ctf_zero_angstrom']=np.sqrt(10000*wave)
result['defocus_pi_error_nm']=resolution**2/wave/10
result['cs_phase_in_pi']=.5*2.7e7*wave**3/resolution**4
result['delocalization_angstrom']=wave*2.2e4/resolution
result['dqe_ssnr_outputs']=[.1*.5,.1*.08]
result['dqe_count_ratio']=.5/.08
result['translation_amplitude']=[float(np.exp(-2*np.pi**2*k**2)) for k in [1/6,1/3]]
angle_sd=np.deg2rad(1.7); radius=65
result['rotation_tangential_sd_pixels']=radius*angle_sd
result['rotation_edge_tangential_amplitude']=np.exp(-2*np.pi**2*(radius*angle_sd/6)**2)
result['independent_reference_distance_threshold']=16900*(1/50-1/1000)
A=np.array([[1.,0.],[0.,1.],[1.,1.]])
result['rectangular_matrix_singular_values']=np.linalg.svd(A,compute_uv=False).tolist()
assert np.linalg.matrix_rank(A)==2
# A full normal matrix has coupling: simple per-coordinate division is not exact.
y=np.array([1.,2.,3.]); normal=A.T@A
v=np.linalg.solve(normal,A.T@y)
assert np.allclose(v,[1,2])
assert not np.allclose((A.T@y)/np.diag(normal),v)
result['coupled_least_squares_solution']=v.tolist()
H=np.array([.9,.2,-.6]); precision=H@H
mean=precision/(precision+1); variance=precision/(precision+1)**2
result['ctf_merge']={'precision':precision,'unregularized_mse':1/precision,
                     'map_mean':mean,'map_variance':variance,'map_bias_squared':(mean-1)**2,
                     'map_mse':variance+(mean-1)**2}
assert np.isclose(result['ctf_merge']['map_mse'],.4524886877836245)
result['fsc_half_to_full_ssnr']={str(r):2*r/(1-r) for r in [.143,1/7,.2,.5]}
assert np.isclose(result['fsc_half_to_full_ssnr'][str(1/7)],1/3)
# Shared nonstructural component acts like common signal in a half-map correlation.
signal=1.; bias=3.; noise=4.
result['shared_bias_correlation']=(signal+bias)/(signal+bias+noise)
result['unbiased_correlation']=signal/(signal+noise)
assert result['shared_bias_correlation']>result['unbiased_correlation']
path=Path(__file__).with_name('science_checks.json')
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
