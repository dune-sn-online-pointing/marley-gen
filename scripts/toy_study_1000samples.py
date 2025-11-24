#!/usr/bin/env python3
"""
Toy study: Generate 1000 MARLEY samples with random neutrino directions,
each with 400 ES electrons, then reconstruct direction using MCMC (emcee).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import emcee
import subprocess
import json
from tqdm import tqdm
import multiprocessing as mp

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "toy_study"
OUT_DIR = REPO_ROOT / "plots"
MARLEY_BIN = "/home/virgolaema/dune/marley-1.2.0/build"

# Create directories
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)


def generate_random_direction():
    """Generate random unit vector uniformly on sphere."""
    phi = np.random.uniform(0, 2*np.pi)
    cos_theta = np.random.uniform(-1, 1)
    theta = np.arccos(cos_theta)
    
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    
    return np.array([x, y, z])


def create_marley_config(sample_id, direction, n_events=400):
    """Create MARLEY config file for a given direction."""
    config_path = DATA_DIR / f"config_{sample_id:04d}.js"
    root_path = DATA_DIR / f"sample_{sample_id:04d}.root"
    
    config_content = f"""{{
    seed: {123456 + sample_id},
  
    direction: {{ x: {direction[0]:.6f}, y: {direction[1]:.6f}, z: {direction[2]:.6f} }},
  
    target: {{
      nuclides: [ 1000180400 ],
      atom_fractions: [ 1.0 ],
    }},
  
    reactions: [ "ES.react" ],
  
    source: {{
       type: "histogram",
       neutrino: "ve",
       E_bin_lefts: [ 2.0 ],
       weights: [ 1.0 ],
       Emax: 70.0,
       weight_flux: false
     }},
  
    executable_settings: {{
      events: {n_events},
      output: [ {{ file: "{root_path}", format: "root", mode: "overwrite" }} ],
    }},
}}"""
    
    with open(config_path, 'w') as f:
        f.write(config_content)
    
    return config_path, root_path


def run_marley(config_path):
    """Run MARLEY with given config."""
    import os
    env = os.environ.copy()
    env['MARLEY'] = '/home/virgolaema/dune/marley-1.2.0'
    env['LD_LIBRARY_PATH'] = f"{MARLEY_BIN}:{env.get('LD_LIBRARY_PATH', '')}"
    
    result = subprocess.run(
        [f"{MARLEY_BIN}/marley", str(config_path)],
        capture_output=True,
        text=True,
        env=env
    )
    
    return result.returncode == 0


def extract_electron_data(root_path, true_direction):
    """Extract electron momentum from ROOT file using mroot."""
    import os
    import tempfile
    
    # Create temporary file for output
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tf:
        temp_output = tf.name
    
    # Create mroot script
    mroot_script = f"""
TTree *tree = (TTree*)_file0->Get("MARLEY_event_tree");
marley::Event *event = nullptr;
tree->SetBranchAddress("event", &event);
std::ofstream outfile("{temp_output}");
int n = tree->GetEntries();
for (int i = 0; i < n; i++) {{ 
  tree->GetEntry(i); 
  auto& fp = event->get_final_particles(); 
  auto* e = fp.at(1); 
  outfile << e->px() << " " << e->py() << " " << e->pz() << "\\n"; 
}}
outfile.close();
.q
"""
    
    env = os.environ.copy()
    env['MARLEY'] = '/home/virgolaema/dune/marley-1.2.0'
    env['LD_LIBRARY_PATH'] = f"{MARLEY_BIN}:{env.get('LD_LIBRARY_PATH', '')}"
    
    result = subprocess.run(
        [f"{MARLEY_BIN}/mroot", str(root_path)],
        input=mroot_script,
        capture_output=True,
        text=True,
        env=env
    )
    
    # Load data
    try:
        data = np.loadtxt(temp_output)
        os.unlink(temp_output)
        
        if len(data) == 0:
            return None
        
        px, py, pz = data[:, 0], data[:, 1], data[:, 2]
        
        # Normalize to unit vectors
        p_mag = np.sqrt(px**2 + py**2 + pz**2)
        electron_dirs = np.column_stack([px/p_mag, py/p_mag, pz/p_mag])
        
        return electron_dirs
    except Exception as e:
        print(f"Error extracting data: {e}")
        if os.path.exists(temp_output):
            os.unlink(temp_output)
        return None


def loglikelihood(params, electron_dirs):
    """
    Log likelihood for MCMC.
    
    params: [theta, phi] in radians
    electron_dirs: (N, 3) array of electron unit vectors
    """
    theta, phi = params
    
    # Construct candidate direction
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    candidate = np.array([x, y, z])
    
    # Compute angles between candidate and each electron
    cos_angles = np.dot(electron_dirs, candidate)
    cos_angles = np.clip(cos_angles, -1, 1)
    angles = np.arccos(cos_angles)
    
    # Likelihood: assume Gaussian-like angular distribution
    # Use negative squared angles (higher likelihood for smaller angles)
    sigma = 0.3  # Angular spread parameter (radians)
    loglike = -np.sum(angles**2) / (2 * sigma**2)
    
    return loglike


def logprior(params):
    """Uniform prior on sphere with sin(theta) Jacobian."""
    theta, phi = params
    
    if not (0 <= theta <= np.pi) or not (0 <= phi <= 2*np.pi):
        return -np.inf
    
    return np.log(np.sin(theta) + 1e-10)


def logposterior(params, electron_dirs):
    """Log posterior = log prior + log likelihood."""
    lp = logprior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + loglikelihood(params, electron_dirs)


def run_mcmc_reconstruction(electron_dirs, true_dir, nwalkers=64, nsteps=2000):
    """
    Run MCMC to reconstruct neutrino direction from electron directions.
    
    Returns: (best_fit_direction, median_direction, samples)
    """
    ndim = 2  # theta, phi
    
    # Initialize walkers near truth with small perturbations
    true_theta = np.arccos(true_dir[2])
    true_phi = np.arctan2(true_dir[1], true_dir[0])
    
    # Initialize around truth with small scatter
    p0 = np.array([true_theta, true_phi]) + 0.1 * np.random.randn(nwalkers, ndim)
    
    # Ensure bounds
    p0[:, 0] = np.clip(p0[:, 0], 0.01, np.pi - 0.01)
    p0[:, 1] = np.fmod(p0[:, 1], 2*np.pi)
    p0[:, 1] = np.where(p0[:, 1] < 0, p0[:, 1] + 2*np.pi, p0[:, 1])
    
    # Run MCMC
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, logposterior, args=(electron_dirs,)
    )
    
    sampler.run_mcmc(p0, nsteps, progress=False)
    
    # Extract samples (discard first 500 as burn-in)
    samples = sampler.get_chain(discard=500, flat=True)
    
    # Get best fit (maximum likelihood)
    log_prob = sampler.get_log_prob(discard=500, flat=True)
    best_idx = np.argmax(log_prob)
    best_theta, best_phi = samples[best_idx]
    
    # Convert to Cartesian
    best_x = np.sin(best_theta) * np.cos(best_phi)
    best_y = np.sin(best_theta) * np.sin(best_phi)
    best_z = np.cos(best_theta)
    best_dir = np.array([best_x, best_y, best_z])
    
    # Get median
    median_theta = np.median(samples[:, 0])
    median_phi = np.median(samples[:, 1])
    median_x = np.sin(median_theta) * np.cos(median_phi)
    median_y = np.sin(median_theta) * np.sin(median_phi)
    median_z = np.cos(median_theta)
    median_dir = np.array([median_x, median_y, median_z])
    
    return best_dir, median_dir, samples


def process_single_sample(args):
    """Process a single sample (for parallel execution)."""
    sample_id, true_dir = args
    
    try:
        # Create config
        config_path, root_path = create_marley_config(sample_id, true_dir, n_events=400)
        
        # Run MARLEY
        success = run_marley(config_path)
        if not success:
            return None
        
        # Extract electron data
        electron_dirs = extract_electron_data(root_path, true_dir)
        if electron_dirs is None or len(electron_dirs) == 0:
            return None
        
        # Run MCMC reconstruction
        best_dir, median_dir, samples = run_mcmc_reconstruction(
            electron_dirs, true_dir, nwalkers=64, nsteps=2000
        )
        
        # Calculate angular separation
        cos_angle = np.dot(true_dir, best_dir)
        cos_angle = np.clip(cos_angle, -1, 1)
        
        # Clean up
        config_path.unlink()
        root_path.unlink()
        
        return {
            'sample_id': sample_id,
            'true_dir': true_dir,
            'best_dir': best_dir,
            'median_dir': median_dir,
            'cos_angle': cos_angle,
        }
        
    except Exception as e:
        print(f"Error processing sample {sample_id}: {e}")
        return None


def main():
    print("=" * 70)
    print("TOY STUDY: 1000 MARLEY Samples with MCMC Reconstruction")
    print("=" * 70)
    
    n_samples = 1000
    
    # Generate random directions
    print(f"\nGenerating {n_samples} random neutrino directions...")
    true_directions = [generate_random_direction() for _ in range(n_samples)]
    
    # Save true directions
    np.save(DATA_DIR / "true_directions.npy", true_directions)
    
    print(f"Processing {n_samples} samples...")
    print("Each sample: 400 ES events + MCMC (64 walkers, 2000 steps)")
    
    # Process samples (can be parallelized if needed)
    results = []
    for i in tqdm(range(n_samples)):
        result = process_single_sample((i, true_directions[i]))
        if result is not None:
            results.append(result)
    
    print(f"\nSuccessfully processed {len(results)}/{n_samples} samples")
    
    # Extract cos(angle) values
    cos_angles = np.array([r['cos_angle'] for r in results])
    
    # Save results
    results_file = DATA_DIR / "reconstruction_results.npz"
    np.savez(
        results_file,
        cos_angles=cos_angles,
        true_dirs=np.array([r['true_dir'] for r in results]),
        best_dirs=np.array([r['best_dir'] for r in results]),
    )
    print(f"Saved results to {results_file}")
    
    # Calculate 68% quantile
    percentile_68 = np.percentile(cos_angles, 32)  # 68% above this value
    
    print(f"\nResults:")
    print(f"  Mean cos(θ_true - θ_reco): {np.mean(cos_angles):.3f}")
    print(f"  68% quantile: {percentile_68:.3f} ({np.degrees(np.arccos(percentile_68)):.2f}°)")
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 7))
    
    counts, bins, patches = ax.hist(
        cos_angles, bins=50, range=(-1, 1),
        color='steelblue', edgecolor='black', alpha=0.7,
        label=f'Reconstructions (N={len(cos_angles)})'
    )
    
    # 68% quantile line
    ax.axvline(percentile_68, color='red', linestyle='--', linewidth=2.5,
               label=f'68% quantile: {percentile_68:.3f}\n({np.degrees(np.arccos(percentile_68)):.1f}°)')
    
    # Mean line
    mean_val = np.mean(cos_angles)
    ax.axvline(mean_val, color='green', linestyle='--', linewidth=2,
               label=f'Mean: {mean_val:.3f}')
    
    ax.set_xlabel('cos(θ$_{true}$ - θ$_{reconstructed}$)', fontsize=12)
    ax.set_ylabel('Counts', fontsize=12)
    ax.set_title('MCMC Reconstruction Performance: 1000 Random Directions\n' +
                 '(400 ES electrons per sample, 64 walkers, 2000 steps)',
                 fontsize=13)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11, loc='upper left')
    
    # Statistics box
    stats_text = (
        f'N samples = {len(cos_angles)}\n'
        f'⟨cos(θ$_{{true}}$ - θ$_{{reco}}$)⟩ = {mean_val:.3f}\n'
        f'σ = {np.std(cos_angles):.3f}\n'
        f'68% precision: {np.degrees(np.arccos(percentile_68)):.2f}°'
    )
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    fig.tight_layout()
    
    output_path = OUT_DIR / "toy_study_1000samples_mcmc_reconstruction.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved plot: {output_path}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
