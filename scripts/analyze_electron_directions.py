#!/usr/bin/python3
"""
Analyze electron directions relative to neutrino direction for each SNB burst.
Computes the average electron direction per burst and the angle to the neutrino direction.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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

def main():
    data_dir = Path("/afs/cern.ch/work/e/evilla/private/dune/marley-gen/data/snb_1000")
    
    print("Computing angles between neutrino and average electron directions...")
    cos_angles = compute_angle_cosines(data_dir)
    
    print(f"\nAnalyzed {len(cos_angles)} bursts")
    print(f"cos(angle) statistics:")
    print(f"  Mean: {cos_angles.mean():.4f}")
    print(f"  Std: {cos_angles.std():.4f}")
    print(f"  Min: {cos_angles.min():.4f}")
    print(f"  Max: {cos_angles.max():.4f}")
    
    # Compute 68% quantile
    quantile_68 = np.quantile(cos_angles, 0.68)
    print(f"\n68% quantile: {quantile_68:.4f}")
    print(f"  (68% of bursts have cos(angle) ≤ {quantile_68:.4f})")
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Histogram
    counts, bins, patches = ax.hist(cos_angles, bins=50, alpha=0.7, 
                                     edgecolor='black', linewidth=0.5,
                                     label=f'N = {len(cos_angles)} bursts')
    
    # Add vertical line for 68% quantile
    ax.axvline(quantile_68, color='red', linestyle='--', linewidth=2,
               label=f'68% quantile: {quantile_68:.4f}')
    
    # Labels and formatting
    ax.set_xlabel('cos(θ) [neutrino dir · avg electron dir]', fontsize=12)
    ax.set_ylabel('Number of bursts', fontsize=12)
    ax.set_title('Angular Distribution: Neutrino Direction vs Average Electron Direction', 
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Save figure
    output_file = data_dir / "electron_direction_analysis.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved to: {output_file}")
    
    # Save data
    output_data = data_dir / "electron_direction_analysis.npz"
    np.savez(output_data, 
             cos_angles=cos_angles,
             quantile_68=quantile_68,
             mean=cos_angles.mean(),
             std=cos_angles.std())
    print(f"✓ Data saved to: {output_data}")
    
    plt.show()

if __name__ == "__main__":
    main()
