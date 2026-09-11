"""Reproduce the original SPA teaching figures and their numeric examples."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'book/images/spa'
EVIDENCE = ROOT / 'notes/spa_revision_20260911'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 12, 'axes.spines.top': False,
                         'axes.spines.right': False, 'figure.dpi': 160})
    colors = ['#246a9e', '#d07a27']
    delta = np.linspace(-8, 8, 301)
    variances = [1., 4.]
    probs = [1 / (1 + np.exp(-2 / (2 * s))) for s in variances]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    for j, (s, p) in enumerate(zip(variances, probs)):
        axs[0].bar(np.array([0, 1]) + (j - .5) * .32, [p, 1-p], width=.30,
                   color=colors[j], label=rf'$\sigma^2={s:g}$')
        for k,v in enumerate([p,1-p]):
            axs[0].text(k+(j-.5)*.32,v+.025,f'{v:.3f}',ha='center',fontsize=10)
        axs[1].plot(delta, 1/(1+np.exp(-delta/(2*s))), color=colors[j],
                    label=rf'$\sigma^2={s:g}$', lw=2.5)
    axs[0].set(xticks=[0,1], xticklabels=['Candidate A', 'Candidate B'], ylim=(0,1),
               ylabel='Posterior probability', title=r'$E_A=10,\ E_B=12$; equal priors')
    axs[0].legend(frameon=False)
    axs[1].axhline(.5,color='#999999',ls=':',lw=1)
    axs[1].set(xlabel=r'Residual difference $E_B-E_A$', ylabel='Posterior for candidate A',
               ylim=(0,1),title='More noise, less decisive weights')
    axs[1].legend(frameon=False)
    fig.savefig(OUT/'em_uncertainty.png'); plt.close(fig)

    rng=np.random.default_rng(20260911)
    covariance=np.array([[5.,4.],[4.,5.]])
    samples=rng.multivariate_normal([0,0],covariance,size=180)
    u=np.array([1.,1.])/np.sqrt(2)
    projections=samples@u[:,None]@u[None,:]
    energy=np.array([9.,4.,.2,0.]); ranks=np.arange(1,5)
    bias=np.array([energy[r:].sum() for r in ranks]); noise=ranks.astype(float)
    fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    axs[0].scatter(*samples.T,s=12,color=colors[0],alpha=.45,label='Noisy images')
    axs[0].plot([-7,7],[-7,7],color=colors[1],lw=2,label='First principal axis')
    for a,b in zip(samples[:15],projections[:15]):
        axs[0].plot([a[0],b[0]],[a[1],b[1]],color='#a8a8a8',lw=.8)
    axs[0].set(xlabel='Pixel 1 (centered)',ylabel='Pixel 2 (centered)',
               xlim=(-7,7),ylim=(-7,7),title='Shared variation in two pixels',aspect='equal')
    axs[0].legend(frameon=False,fontsize=10)
    axs[1].plot(ranks,bias,'o-',color=colors[0],label='Discarded signal')
    axs[1].plot(ranks,noise,'s-',color=colors[1],label='Retained noise')
    axs[1].plot(ranks,bias+noise,'D-',color='#416d48',lw=2,label='Expected squared error')
    axs[1].set(xticks=ranks,xlabel='Retained rank',ylabel='Energy / squared error',
               title=r'Known fixed basis; $\sigma^2=1$',ylim=(0,6))
    axs[1].legend(frameon=False,fontsize=10)
    fig.savefig(OUT/'pca_geometry.png');plt.close(fig)
    frequency=np.linspace(0,.25,401)
    discrete=np.abs(np.cos(np.pi*5*frequency))
    continuous=np.exp(-2*np.pi**2*2**2*frequency**2)
    fig,ax=plt.subplots(figsize=(8,3.8),layout='constrained')
    ax.plot(frequency,discrete,color=colors[0],lw=2.5,label=r'Two equal states; $d=5$ Å')
    ax.plot(frequency,continuous,color=colors[1],lw=2.5,label=r'Gaussian displacement; $\sigma_d=2$ Å')
    ax.scatter([.1,.2],[0,1],color=colors[0],zorder=3)
    ax.set(xlabel='Spatial frequency along displacement (Å⁻¹)',ylabel='Retained amplitude',
           title='A discrete mixture oscillates; Gaussian motion damps',ylim=(-.03,1.1),xlim=(0,.25))
    ax.legend(frameon=False,fontsize=10,loc='center right')
    fig.savefig(OUT/'heterogeneity_attenuation.png');plt.close(fig)
    np.savetxt(EVIDENCE/'heterogeneity_attenuation.csv',np.column_stack([frequency,discrete,continuous]),
               delimiter=',',header='frequency_inverse_angstrom,two_states_amplitude,gaussian_amplitude',comments='')
    x=np.array([0.,1.,3.,4.]); means=np.array([0.,4.])
    weights=np.exp(-.5*(x[:,None]-means[None,:])**2)
    weights/=weights.sum(axis=1,keepdims=True)
    updated=(weights*x[:,None]).sum(axis=0)/weights.sum(axis=0)
    assert np.allclose(updated,[.518656910,.0+3.481343090],atol=1e-8)
    assert np.allclose(bias+noise,[5.2,2.2,3,4])
    assert np.allclose(np.linalg.eigvalsh(covariance),[1.,9.])
    report={'seed':20260911,'posterior_pose_A':probs,'em_responsibilities':weights.tolist(),
            'em_updated_means':updated.tolist(),'pca_population_covariance':covariance.tolist(),
            'fixed_basis_signal_energies':energy.tolist(),'rank_risk':(bias+noise).tolist(),
            'heterogeneity_example':{'state_separation_angstrom':5,'gaussian_sd_angstrom':2},
            'map_gain':.1/(.1**2+.04),'shrinkage_mse':(.2-1)**2+.2**2*4,
            'outputs':['book/images/spa/em_uncertainty.png','book/images/spa/pca_geometry.png','book/images/spa/heterogeneity_attenuation.png']}
    (EVIDENCE/'numerical_examples.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savetxt(EVIDENCE/'pca_illustration_points.csv',samples,delimiter=',',header='pixel1,pixel2',comments='')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
