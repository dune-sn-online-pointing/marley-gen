#!/usr/bin/python3
"""
MCMC analysis of electron direction distribution using emcee.
Fits the cos(angle) distribution between neutrino and average electron directions.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import emcee
import corner

def compute_angle_cosines(data_dir):
    """
    For each burst:
    - Compute average electron direction
    - Compute angle between average electron direction and neutrino direction
    - Return cos(angle) for all bursts
    """
    job_files = sorted(Path(data_dir).glob("job_*_results.npz"))
    
    cos_angles = []
    
    for job_file in job_files:
        data = np.load(job_file, allow_pickle=True)
        directions = data['directions']  # (n_bursts, 3)
        electrons = data['electrons']     # (n_bursts, 400, 3)
        
        for i, (nu_dir, elec_momenta) in enumerate(zip(directions, electrons)):
            if elec_momenta is None or len(elec_momenta) == 0:
                continue
            
            # Compute average electron momentum direction
            avg_elec_momentum = elec_momenta.mean(axis=0)  # (3,)
            
            # Normalize both directions
            nu_dir_norm = nu_dir / np.linalg.norm(nu_dir)
            avg_elec_dir_norm = avg_elec_momentum / np.linalg.norm(avg_elec_momentum)
            
            # Compute cos(angle) = dot product of normalized vectors
            cos_angle = np.dot(nu_dir_norm, avg_elec_dir_norm)
            cos_angles.append(cos_angle)
    
    return np.array(cos_angles)

def log_likelihood(theta, cos_angles):
    """
    Log likelihood for a truncated normal distribution on cos(angle).
    Parameters: mu, sigma
    """
    mu, log_sigma = theta
    sigma = np.exp(log_sigma)
    
    # Truncated normal: valid range is [0, 1] for cos(angle)
    # But our data is very close to 1, so we model it as normal
    model = -0.5 * np.sum(((cos_angles - mu) / sigma) ** 2)
    model -= len(cos_angles) * np.log(sigma)
    
    return model

def log_prior(theta):
    """
    Log prior for parameters.
    mu: mean of cos(angle), should be close to 1
    log_sigma: log of standard deviation
    """
    mu, log_sigma = theta
    
    # Prior on mu: should be between 0.99 and 1.0
    if not (0.99 < mu < 1.001):
        return -np.inf
    
    # Prior on log_sigma: reasonable range for std
    if not (-10 < log_sigma < -3):
        return -np.inf
    
    return 0.0

def log_probability(theta, cos_angles):
    """Log posterior probability."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, cos_angles)

def run_mcmc(cos_angles, n_walkers=128, n_steps=200):
    """
    Run MCMC sampling using emcee.
    """
    # Initialize walkers around the data mean and std
    mu_init = cos_angles.mean()
    sigma_init = cos_angles.std()
    
    print(f"Initial estimates: mu={mu_init:.6f}, sigma={sigma_init:.6f}")
    
    # Parameter space: [mu, log_sigma]
    ndim = 2
    
    # Initialize walkers in a small ball around initial guess
    pos = np.zeros((n_walkers, ndim))
    pos[:, 0] = mu_init + 1e-5 * np.random.randn(n_walkers)
    pos[:, 1] = np.log(sigma_init) + 0.1 * np.random.randn(n_walkers)
    
    # Set up the sampler
    sampler = emcee.EnsembleSampler(
        n_walkers, ndim, log_probability, args=[cos_angles]
    )
    
    # Run MCMC
    print(f"\nRunning MCMC with {n_walkers} walkers for {n_steps} steps...")
    sampler.run_mcmc(pos, n_steps, progress=True)
    
    return sampler

def analyze_chains(sampler, burnin=50):
    """
    Analyze MCMC chains and compute statistics.
    """
    # Get the chain
    samples = sampler.get_chain()
    
    print(f"\nChain shape: {samples.shape}")
    print(f"Acceptance fraction: {np.mean(sampler.acceptance_fraction):.3f}")
    
    # Discard burn-in
    flat_samples = sampler.get_chain(discard=burnin, thin=1, flat=True)
    
    # Transform log_sigma back to sigma
    flat_samples[:, 1] = np.exp(flat_samples[:, 1])
    
    # Compute statistics
    mu_mean = np.mean(flat_samples[:, 0])
    mu_std = np.std(flat_samples[:, 0])
    sigma_mean = np.mean(flat_samples[:, 1])
    sigma_std = np.std(flat_samples[:, 1])
    
    # Compute quantiles
    mu_16, mu_50, mu_84 = np.percentile(flat_samples[:, 0], [16, 50, 84])
    sigma_16, sigma_50, sigma_84 = np.percentile(flat_samples[:, 1], [16, 50, 84])
    
    print(f"\nMCMC Results (after {burnin} step burn-in):")
    print(f"  μ = {mu_50:.6f} +{mu_84-mu_50:.6f} -{mu_50-mu_16:.6f}")
    print(f"  σ = {sigma_50:.6f} +{sigma_84-sigma_50:.6f} -{sigma_50-sigma_16:.6f}")
    print(f"\nMean ± Std:")
    print(f"  μ = {mu_mean:.6f} ± {mu_std:.6f}")
    print(f"  σ = {sigma_mean:.6f} ± {sigma_std:.6f}")
    
    return flat_samples, samples

def plot_results(cos_angles, flat_samples, samples, output_dir, burnin=50):
    """
    Create diagnostic plots.
    """
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Trace plots
    ax1 = plt.subplot(3, 3, 1)
    ax1.plot(samples[:, :, 0], 'k', alpha=0.1)
    ax1.axvline(burnin, color='r', linestyle='--', label='Burn-in')
    ax1.set_xlabel('Step')
    ax1.set_ylabel('μ (cos θ)')
    ax1.set_title('Trace: μ')
    ax1.legend()
    
    ax2 = plt.subplot(3, 3, 2)
    ax2.plot(np.exp(samples[:, :, 1]), 'k', alpha=0.1)
    ax2.axvline(burnin, color='r', linestyle='--', label='Burn-in')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('σ')
    ax2.set_title('Trace: σ')
    ax2.legend()
    
    # 2. Corner plot
    ax3 = plt.subplot(3, 3, 3)
    mu_samples = flat_samples[:, 0]
    sigma_samples = flat_samples[:, 1]
    
    # 2D histogram
    h, xedges, yedges = np.histogram2d(mu_samples, sigma_samples, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    im = ax3.imshow(h.T, origin='lower', extent=extent, aspect='auto', cmap='Blues')
    ax3.set_xlabel('μ (cos θ)')
    ax3.set_ylabel('σ')
    ax3.set_title('Joint Posterior')
    plt.colorbar(im, ax=ax3)
    
    # 3. Posterior distributions
    ax4 = plt.subplot(3, 3, 4)
    ax4.hist(mu_samples, bins=50, alpha=0.7, edgecolor='black', density=True)
    mu_16, mu_50, mu_84 = np.percentile(mu_samples, [16, 50, 84])
    ax4.axvline(mu_50, color='r', linestyle='-', linewidth=2, label=f'Median: {mu_50:.6f}')
    ax4.axvline(mu_16, color='r', linestyle='--', alpha=0.5)
    ax4.axvline(mu_84, color='r', linestyle='--', alpha=0.5)
    ax4.set_xlabel('μ (cos θ)')
    ax4.set_ylabel('Posterior density')
    ax4.set_title('Posterior: μ')
    ax4.legend()
    
    ax5 = plt.subplot(3, 3, 5)
    ax5.hist(sigma_samples, bins=50, alpha=0.7, edgecolor='black', density=True)
    sigma_16, sigma_50, sigma_84 = np.percentile(sigma_samples, [16, 50, 84])
    ax5.axvline(sigma_50, color='r', linestyle='-', linewidth=2, label=f'Median: {sigma_50:.6f}')
    ax5.axvline(sigma_16, color='r', linestyle='--', alpha=0.5)
    ax5.axvline(sigma_84, color='r', linestyle='--', alpha=0.5)
    ax5.set_xlabel('σ')
    ax5.set_ylabel('Posterior density')
    ax5.set_title('Posterior: σ')
    ax5.legend()
    
    # 4. Data vs Model
    ax6 = plt.subplot(3, 3, 6)
    counts, bins, _ = ax6.hist(cos_angles, bins=50, alpha=0.7, edgecolor='black', 
                                label=f'Data (N={len(cos_angles)})', density=True)
    
    # Overlay model samples
    x = np.linspace(cos_angles.min(), cos_angles.max(), 200)
    for mu, sigma in flat_samples[np.random.randint(len(flat_samples), size=100)]:
        model = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
        ax6.plot(x, model, 'r-', alpha=0.02)
    
    # Best fit
    mu_best = mu_50
    sigma_best = sigma_50
    model_best = (1 / (sigma_best * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu_best) / sigma_best) ** 2)
    ax6.plot(x, model_best, 'b-', linewidth=2, label=f'Best fit: μ={mu_best:.5f}, σ={sigma_best:.5f}')
    
    ax6.set_xlabel('cos(θ) [neutrino dir · avg electron dir]')
    ax6.set_ylabel('Density')
    ax6.set_title('Data vs MCMC Model')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    # 5. Quantiles comparison
    ax7 = plt.subplot(3, 3, 7)
    quantiles = [0.16, 0.50, 0.68, 0.84, 0.95]
    data_quantiles = np.percentile(cos_angles, np.array(quantiles) * 100)
    
    # Model quantiles from posterior samples
    model_quantiles = []
    for q in quantiles:
        q_samples = []
        for mu, sigma in flat_samples[np.random.randint(len(flat_samples), size=1000)]:
            # For normal distribution, quantile is mu + sigma * z
            from scipy.stats import norm
            z = norm.ppf(q)
            q_samples.append(mu + sigma * z)
        model_quantiles.append(np.mean(q_samples))
    
    ax7.plot(quantiles, data_quantiles, 'bo-', label='Data', markersize=8)
    ax7.plot(quantiles, model_quantiles, 'rs-', label='MCMC Model', markersize=8)
    ax7.set_xlabel('Quantile')
    ax7.set_ylabel('cos(θ)')
    ax7.set_title('Quantile Comparison')
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    
    # 6. Autocorrelation
    ax8 = plt.subplot(3, 3, 8)
    try:
        tau = sampler.get_autocorr_time()
        ax8.bar([0, 1], tau, color=['blue', 'red'], alpha=0.7)
        ax8.set_xticks([0, 1])
        ax8.set_xticklabels(['μ', 'σ'])
        ax8.set_ylabel('Autocorrelation time')
        ax8.set_title(f'Autocorrelation Time\nμ: {tau[0]:.1f}, σ: {tau[1]:.1f}')
        ax8.grid(True, alpha=0.3, axis='y')
    except Exception as e:
        ax8.text(0.5, 0.5, f'Autocorrelation\ncalculation failed:\n{str(e)}', 
                ha='center', va='center', transform=ax8.transAxes)
        ax8.set_title('Autocorrelation Time')
    
    # 7. Residuals
    ax9 = plt.subplot(3, 3, 9)
    mu_best = mu_50
    sigma_best = sigma_50
    expected_quantiles = np.linspace(0, 1, len(cos_angles))
    observed_quantiles = np.searchsorted(np.sort(cos_angles), cos_angles) / len(cos_angles)
    
    # Q-Q plot
    from scipy.stats import norm
    theoretical = norm.ppf(expected_quantiles[1:-1], loc=mu_best, scale=sigma_best)
    observed = np.sort(cos_angles)[1:-1]
    ax9.scatter(theoretical, observed, alpha=0.5, s=10)
    ax9.plot([observed.min(), observed.max()], [observed.min(), observed.max()], 
            'r--', linewidth=2, label='Perfect fit')
    ax9.set_xlabel('Theoretical quantiles')
    ax9.set_ylabel('Observed quantiles')
    ax9.set_title('Q-Q Plot')
    ax9.legend()
    ax9.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = output_dir / "mcmc_electron_direction_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ MCMC diagnostic plot saved to: {output_file}")
    
    # Create corner plot separately
    fig2 = corner.corner(
        flat_samples,
        labels=['μ (cos θ)', 'σ'],
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_fmt='.6f'
    )
    
    corner_file = output_dir / "mcmc_corner_plot.png"
    fig2.savefig(corner_file, dpi=300, bbox_inches='tight')
    print(f"✓ Corner plot saved to: {corner_file}")

def main():
    data_dir = Path("/afs/cern.ch/work/e/evilla/private/dune/marley-gen/data/snb_1000")
    
    print("Loading data and computing angles...")
    cos_angles = compute_angle_cosines(data_dir)
    
    print(f"\nAnalyzed {len(cos_angles)} bursts")
    print(f"cos(angle) statistics:")
    print(f"  Mean: {cos_angles.mean():.6f}")
    print(f"  Std: {cos_angles.std():.6f}")
    print(f"  Min: {cos_angles.min():.6f}")
    print(f"  Max: {cos_angles.max():.6f}")
    
    # Run MCMC
    n_walkers = 128
    n_steps = 200
    burnin = 50
    
    sampler = run_mcmc(cos_angles, n_walkers=n_walkers, n_steps=n_steps)
    
    # Analyze results
    flat_samples, samples = analyze_chains(sampler, burnin=burnin)
    
    # Create plots
    plot_results(cos_angles, flat_samples, samples, data_dir, burnin=burnin)
    
    # Save results
    output_data = data_dir / "mcmc_electron_direction_results.npz"
    np.savez(
        output_data,
        cos_angles=cos_angles,
        flat_samples=flat_samples,
        samples=samples,
        n_walkers=n_walkers,
        n_steps=n_steps,
        burnin=burnin
    )
    print(f"✓ MCMC results saved to: {output_data}")

if __name__ == "__main__":
    main()
