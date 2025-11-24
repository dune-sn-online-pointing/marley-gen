#!/usr/bin/env python3
"""Plot the angular distribution of electrons from MARLEY ES events.

Handles both old (flat branches) and new (TTree with objects) ROOT formats.
"""
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt
import uproot
import awkward as ak

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "plots"


def load_electron_angles(root_file):
    """Load electron momentum components and calculate angles.
    
    Returns:
        dict with:
        - px, py, pz: momentum components (MeV/c)
        - theta: polar angle (radians)
        - phi: azimuthal angle (radians)
        - cos_theta: cosine of polar angle
    """
    root_path = Path(root_file)
    if not root_path.is_absolute():
        root_path = (REPO_ROOT / root_path).resolve()

    if not root_path.exists():
        raise FileNotFoundError(f"MARLEY file not found: {root_path}")

    with uproot.open(root_path) as f:
        # Try old format first (flat branches)
        if 'mst' in f:
            tree = f["mst"]
            px = tree["pxr"].array(library="np")
            py = tree["pyr"].array(library="np")
            pz = tree["pzr"].array(library="np")
        
        # New format (object branches)
        elif 'MARLEY_event_tree' in f:
            tree = f["MARLEY_event_tree"]
            events = tree['event'].array(library="ak")
            
            # Extract final state particles (index 1 is typically the electron for ES)
            # ES: nu + e -> nu' + e', so final_particles[0]=nu', final_particles[1]=e'
            final_particles = events.final_particles_
            
            # Get electron (should be second final particle, index 1)
            electrons = final_particles[:, 1]
            
            px = ak.to_numpy(electrons.px_)
            py = ak.to_numpy(electrons.py_)
            pz = ak.to_numpy(electrons.pz_)
        else:
            raise ValueError(f"Unknown ROOT file format. Available keys: {list(f.keys())}")
        
        # Calculate angles
        p_mag = np.sqrt(px**2 + py**2 + pz**2)
        
        # Polar angle (theta): angle from +z axis
        cos_theta = np.divide(pz, p_mag, out=np.ones_like(pz), where=p_mag > 0)
        theta = np.arccos(np.clip(cos_theta, -1, 1))  # radians
        
        # Azimuthal angle (phi): angle in x-y plane from +x axis
        phi = np.arctan2(py, px)  # radians, range [-pi, pi]
        
        return {
            'px': px,
            'py': py,
            'pz': pz,
            'p_mag': p_mag,
            'theta': theta,
            'cos_theta': cos_theta,
            'phi': phi,
        }


def plot_angle_distributions(angles, output_dir, tag="es_angles"):
    """Create comprehensive angular distribution plots."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a figure with multiple subplots
    fig = plt.figure(figsize=(14, 10))
    
    # 1. Polar angle (theta) distribution
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.hist(np.degrees(angles['theta']), bins=50, range=(0, 180),
             color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Polar Angle θ (degrees)')
    ax1.set_ylabel('Counts')
    ax1.set_title('Electron Polar Angle Distribution')
    ax1.grid(alpha=0.3)
    ax1.axvline(90, color='red', linestyle='--', alpha=0.5, label='90°')
    ax1.legend()
    
    # 2. cos(theta) distribution
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.hist(angles['cos_theta'], bins=50, range=(-1, 1),
             color='orange', edgecolor='black', alpha=0.7)
    ax2.set_xlabel('cos(θ)')
    ax2.set_ylabel('Counts')
    ax2.set_title('Electron cos(θ) Distribution')
    ax2.grid(alpha=0.3)
    ax2.axvline(0, color='red', linestyle='--', alpha=0.5, label='cos(θ)=0')
    ax2.axvline(1, color='green', linestyle='--', alpha=0.5, label='cos(θ)=1')
    ax2.legend()
    
    # 3. Azimuthal angle (phi) distribution
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.hist(np.degrees(angles['phi']), bins=50, range=(-180, 180),
             color='mediumpurple', edgecolor='black', alpha=0.7)
    ax3.set_xlabel('Azimuthal Angle φ (degrees)')
    ax3.set_ylabel('Counts')
    ax3.set_title('Electron Azimuthal Angle Distribution')
    ax3.grid(alpha=0.3)
    
    # 4. 2D scatter: theta vs phi
    ax4 = fig.add_subplot(2, 3, 4)
    h = ax4.hist2d(np.degrees(angles['phi']), np.degrees(angles['theta']),
                   bins=[40, 40], range=[[-180, 180], [0, 180]],
                   cmap='viridis')
    ax4.set_xlabel('φ (degrees)')
    ax4.set_ylabel('θ (degrees)')
    ax4.set_title('2D Angular Distribution')
    plt.colorbar(h[3], ax=ax4, label='Counts')
    
    # 5. Momentum components
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.hist(angles['pz'], bins=50, alpha=0.6, label='pz', color='blue', edgecolor='black')
    ax5.hist(angles['px'], bins=50, alpha=0.6, label='px', color='red', edgecolor='black')
    ax5.hist(angles['py'], bins=50, alpha=0.6, label='py', color='green', edgecolor='black')
    ax5.set_xlabel('Momentum Component (MeV/c)')
    ax5.set_ylabel('Counts')
    ax5.set_title('Electron Momentum Components')
    ax5.legend()
    ax5.grid(alpha=0.3)
    
    # 6. Polar plot of directions
    ax6 = fig.add_subplot(2, 3, 6, projection='polar')
    # Sample a subset if there are too many events
    n_sample = min(1000, len(angles['theta']))
    indices = np.random.choice(len(angles['theta']), n_sample, replace=False)
    ax6.scatter(angles['phi'][indices], np.degrees(angles['theta'][indices]),
                s=2, alpha=0.5, c='steelblue')
    ax6.set_theta_zero_location('E')
    ax6.set_theta_direction(1)
    ax6.set_ylim(0, 180)
    ax6.set_ylabel('θ (degrees)', labelpad=30)
    ax6.set_title('Polar View of Electron Directions', pad=20)
    
    # Add statistics text
    stats_text = (
        f"Statistics (N={len(angles['theta'])} events):\n"
        f"⟨θ⟩ = {np.degrees(np.mean(angles['theta'])):.2f}°\n"
        f"⟨cos θ⟩ = {np.mean(angles['cos_theta']):.3f}\n"
        f"⟨φ⟩ = {np.degrees(np.mean(angles['phi'])):.2f}°"
    )
    fig.text(0.02, 0.02, stats_text, fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    fig.suptitle('Electron Angular Distributions from ES Events', fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0.05, 1, 0.99])
    
    output_path = output_dir / f"{tag}_comprehensive.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Saved comprehensive angular plot: {output_path}")
    return output_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot electron angular distributions from MARLEY ES events"
    )
    parser.add_argument(
        "--input", "-i", type=Path,
        default=Path("data/flat_ES.root"),
        help="Path to MARLEY ES ROOT file",
    )
    parser.add_argument(
        "--tag", "-t", default="electron_angles",
        help="Tag used in output filename",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    print(f"Loading electron data from {args.input} ...")
    angles = load_electron_angles(args.input)
    
    print(f"Loaded {len(angles['theta'])} events")
    print(f"\nElectron Angular Statistics:")
    print(f"  Mean θ:      {np.degrees(np.mean(angles['theta'])):.2f}°")
    print(f"  Mean cos(θ): {np.mean(angles['cos_theta']):.3f}")
    print(f"  Mean φ:      {np.degrees(np.mean(angles['phi'])):.2f}°")
    print(f"  θ range:     [{np.degrees(angles['theta'].min()):.2f}°, {np.degrees(angles['theta'].max()):.2f}°]")
    
    print("\nCreating angular distribution plots...")
    plot_angle_distributions(angles, OUT_DIR, args.tag)
    
    print("\n✓ Done!")


if __name__ == "__main__":
    main()
