#!/usr/bin/env python3
"""
Create Physics-Based Plots from Flux Data
==========================================
Since we don't have the MARLEY ROOT files, we'll create physics-based plots
using the flux files and known ES/CC kinematics.

For ES: nu + e -> nu' + e'
For CC: nu_e + Ar40 -> e- + K40*
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FLUX_DIR = REPO_ROOT / "interacted-fluxes"

def load_gkvm_flux(mode='ES'):
    """Load GKVM flux data"""
    if mode == 'ES':
        # ES: all flavors contribute
        flavors = ["nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
        files = [FLUX_DIR / f"gvkm_{flav}_e_ar40kt_events.dat" for flav in flavors]
        
        spectra = []
        for f in files:
            if f.exists():
                data = np.genfromtxt(f, delimiter=' ', skip_footer=2)
                spectra.append(data[:, 1])
        
        # Energy is same for all
        data = np.genfromtxt(files[0], delimiter=' ', skip_footer=2)
        energy_gev = data[:, 0]
        energy_mev = energy_gev * 1000
        spectrum = np.sum(spectra, axis=0)
        
    elif mode == 'CC':
        # CC: only nu_e
        f = FLUX_DIR / "gvkm_nue_Ar40_ar40kt_events.dat"
        data = np.genfromtxt(f, delimiter=' ', skip_footer=2)
        energy_gev = data[:, 0]
        energy_mev = energy_gev * 1000
        spectrum = data[:, 1]
    
    return energy_mev, spectrum


def es_kinematics(E_nu, cos_theta_nu):
    """
    Calculate electron recoil energy for ES scattering.
    nu + e -> nu' + e'
    
    Args:
        E_nu: Neutrino energy (MeV)
        cos_theta_nu: Cosine of NEUTRINO scattering angle (can be [-1, 1])
    
    Returns:
        T_e: Electron recoil kinetic energy (MeV)
    
    Formula: T_e = 2·E_ν²·sin²(θ_ν/2) / (m_e + 2·E_ν·sin²(θ_ν/2))
    
    Physical limits:
    - θ_ν = 0° (forward): T_e = 0 (no scattering)
    - θ_ν = 180° (backward): T_e = max (head-on collision)
    """
    m_e = 0.511  # MeV
    
    # Convert cos(θ) to sin²(θ/2) using: sin²(θ/2) = (1 - cos θ)/2
    sin2_half = (1.0 - cos_theta_nu) / 2.0
    
    # ES formula using neutrino scattering angle
    numerator = 2 * E_nu**2 * sin2_half
    denominator = m_e + 2 * E_nu * sin2_half
    T_e = numerator / denominator
    
    return T_e


def cc_kinematics(E_nu):
    """
    For CC interactions, electron takes most of the neutrino energy
    minus the Q-value and some recoil energy.
    
    Approximate: E_e ~ E_nu - Q_value - small recoil
    Q-value for nu_e + Ar40 -> e- + K40* is about 5-10 MeV
    
    Note: CC threshold is ~11 MeV. Events below threshold are filtered out.
    """
    Q_value = 7.0  # MeV (approximate)
    # Add some smearing for nuclear effects
    recoil = np.random.exponential(2.0, size=len(E_nu) if hasattr(E_nu, '__len__') else 1)
    E_e = E_nu - Q_value - recoil
    
    # Don't clamp to minimum - let negative values exist so we can filter them
    return E_e


def generate_es_events_flat(n_events=10000, E_min=1, E_max=100):
    """Generate ES events from flat neutrino spectrum"""
    # Flat spectrum
    E_nu_samples = np.random.uniform(E_min, E_max, n_events)
    
    # Sample neutrino scattering angle from ES differential cross section
    cos_theta_nu_samples = []
    while len(cos_theta_nu_samples) < n_events:
        cos_theta_proposal = np.random.uniform(-1, 1, n_events * 2)
        weight = (1 + cos_theta_proposal)**2 / 4.0
        accept = np.random.uniform(0, 1, len(cos_theta_proposal)) < weight
        cos_theta_nu_samples.extend(cos_theta_proposal[accept])
    
    cos_theta_nu_samples = np.array(cos_theta_nu_samples[:n_events])
    
    # Calculate electron energies
    E_e_samples = es_kinematics(E_nu_samples, cos_theta_nu_samples)
    
    # Calculate electron angles
    E_nu_f = E_nu_samples - E_e_samples
    p_nu_f = E_nu_f
    theta_nu = np.arccos(np.clip(cos_theta_nu_samples, -1, 1))
    p_nu_f_z = p_nu_f * cos_theta_nu_samples
    p_nu_f_x = p_nu_f * np.sin(theta_nu)
    
    p_e_z = E_nu_samples - p_nu_f_z
    p_e_x = -p_nu_f_x
    p_e_mag = np.sqrt(p_e_x**2 + p_e_z**2)
    cos_theta_e_samples = p_e_z / p_e_mag
    
    return {
        'E_nu': E_nu_samples,
        'E_e': E_e_samples,
        'cos_theta_nu': cos_theta_nu_samples,
        'cos_theta_e': cos_theta_e_samples
    }


def generate_cc_events_flat(n_events=10000, E_min=1, E_max=100):
    """Generate CC events from flat neutrino spectrum"""
    # Generate more events than needed to account for threshold filtering
    n_generate = int(n_events * 1.5)
    E_nu_samples = np.random.uniform(E_min, E_max, n_generate)
    
    # Isotropic angular distribution
    cos_theta_samples = np.random.uniform(-1.0, 1.0, size=n_generate)
    
    # Calculate electron energies
    E_e_samples = cc_kinematics(E_nu_samples)
    
    # Filter out events below threshold (E_e < 0)
    valid_mask = E_e_samples > 0
    E_nu_samples = E_nu_samples[valid_mask][:n_events]
    E_e_samples = E_e_samples[valid_mask][:n_events]
    cos_theta_samples = cos_theta_samples[valid_mask][:n_events]
    
    return {
        'E_nu': E_nu_samples,
        'E_e': E_e_samples,
        'cos_theta': cos_theta_samples
    }


def generate_es_events(n_events=10000):
    """Generate ES events from GKVM flux"""
    E_nu_mev, flux = load_gkvm_flux('ES')
    
    # Sample neutrino energies from flux
    flux_norm = flux / np.sum(flux)
    E_nu_samples = np.random.choice(E_nu_mev, size=n_events, p=flux_norm)
    
    # Sample neutrino scattering angle from ES differential cross section
    # dσ/d(cos θ) ∝ (1 + cos θ)^2 for θ = neutrino scattering angle
    # Use rejection sampling to get correct distribution
    cos_theta_nu_samples = []
    while len(cos_theta_nu_samples) < n_events:
        # Propose uniform samples in [-1, 1]
        cos_theta_proposal = np.random.uniform(-1, 1, n_events * 2)
        # Weight by (1 + cos θ)^2, normalized by max value of 4 at cos θ = 1
        weight = (1 + cos_theta_proposal)**2 / 4.0
        accept = np.random.uniform(0, 1, len(cos_theta_proposal)) < weight
        cos_theta_nu_samples.extend(cos_theta_proposal[accept])
    
    cos_theta_nu_samples = np.array(cos_theta_nu_samples[:n_events])
    
    # Calculate electron energies from neutrino scattering angle
    E_e_samples = es_kinematics(E_nu_samples, cos_theta_nu_samples)
    
    # Now calculate ELECTRON angles from momentum conservation
    # p_nu_i = E_nu (along z), p_nu_f at angle theta_nu
    # p_e = p_nu_i - p_nu_f
    E_nu_f = E_nu_samples - E_e_samples  # Energy conservation
    p_nu_f = E_nu_f  # Massless neutrino
    
    # Neutrino momentum components (assuming phi=0)
    theta_nu = np.arccos(np.clip(cos_theta_nu_samples, -1, 1))
    p_nu_f_z = p_nu_f * cos_theta_nu_samples
    p_nu_f_x = p_nu_f * np.sin(theta_nu)
    
    # Electron momentum from conservation
    p_e_z = E_nu_samples - p_nu_f_z
    p_e_x = -p_nu_f_x
    p_e_mag = np.sqrt(p_e_x**2 + p_e_z**2)
    cos_theta_e_samples = p_e_z / p_e_mag
    
    return {
        'E_nu': E_nu_samples,
        'E_e': E_e_samples,
        'cos_theta_nu': cos_theta_nu_samples,  # Neutrino angle
        'cos_theta_e': cos_theta_e_samples     # Electron angle (what we want to plot!)
    }


def generate_cc_events(n_events=10000):
    """Generate CC events from GKVM flux"""
    E_nu_mev, flux = load_gkvm_flux('CC')
    
    # Remove zero flux entries
    mask = flux > 0
    E_nu_mev = E_nu_mev[mask]
    flux = flux[mask]
    
    if len(flux) == 0:
        print("Warning: No non-zero CC flux!")
        return None
    
    # Generate more events than needed to account for threshold filtering
    n_generate = int(n_events * 1.5)
    
    # Sample neutrino energies from flux
    flux_norm = flux / np.sum(flux)
    E_nu_samples = np.random.choice(E_nu_mev, size=n_generate, p=flux_norm)
    
    # For CC, angular distribution is more isotropic
    cos_theta_samples = np.random.uniform(-1.0, 1.0, size=n_generate)
    
    # Calculate electron energies
    E_e_samples = cc_kinematics(E_nu_samples)
    
    # Filter out events below threshold (E_e < 0)
    valid_mask = E_e_samples > 0
    E_nu_samples = E_nu_samples[valid_mask][:n_events]
    E_e_samples = E_e_samples[valid_mask][:n_events]
    cos_theta_samples = cos_theta_samples[valid_mask][:n_events]
    
    return {
        'E_nu': E_nu_samples,
        'E_e': E_e_samples,
        'cos_theta': cos_theta_samples
    }


def plot_electron_energies():
    """Plot electron energy spectra for ES and CC"""
    print("Generating ES events...")
    es_data = generate_es_events(50000)
    
    print("Generating CC events...")
    cc_data = generate_cc_events(50000)
    
    if cc_data is None:
        print("Skipping CC plots - no data")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # ES electron energy - linear
    axes[0, 0].hist(es_data['E_e'], bins=100, range=(0, 70), 
                    alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
    axes[0, 0].set_xlabel('Electron Energy [MeV]', fontsize=12)
    axes[0, 0].set_ylabel('Counts', fontsize=12)
    axes[0, 0].set_title('ES: Electron Energy (GKVM Flux)', fontsize=13, fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)
    
    # ES electron energy - log
    axes[0, 1].hist(es_data['E_e'], bins=100, range=(0, 70),
                    alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
    axes[0, 1].set_xlabel('Electron Energy [MeV]', fontsize=12)
    axes[0, 1].set_ylabel('Counts', fontsize=12)
    axes[0, 1].set_title('ES: Electron Energy (Log Scale)', fontsize=13, fontweight='bold')
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)
    
    # CC electron energy - linear
    axes[1, 0].hist(cc_data['E_e'], bins=100, range=(0, 70),
                    alpha=0.7, color='coral', edgecolor='black', linewidth=0.5)
    axes[1, 0].set_xlabel('Electron Energy [MeV]', fontsize=12)
    axes[1, 0].set_ylabel('Counts', fontsize=12)
    axes[1, 0].set_title('CC: Electron Energy (GKVM Flux)', fontsize=13, fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    
    # CC electron energy - log
    axes[1, 1].hist(cc_data['E_e'], bins=100, range=(0, 70),
                    alpha=0.7, color='coral', edgecolor='black', linewidth=0.5)
    axes[1, 1].set_xlabel('Electron Energy [MeV]', fontsize=12)
    axes[1, 1].set_ylabel('Counts', fontsize=12)
    axes[1, 1].set_title('CC: Electron Energy (Log Scale)', fontsize=13, fontweight='bold')
    axes[1, 1].set_yscale('log')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "electron_energy_es_cc_comparison.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'electron_energy_es_cc_comparison.png'}")
    plt.close()


def plot_cos_theta_distributions():
    """Plot angular distributions"""
    print("Generating ES events for cos(theta)...")
    es_data = generate_es_events(50000)
    
    print("Generating CC events for cos(theta)...")
    cc_data = generate_cc_events(50000)
    
    if cc_data is None:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Separate plots
    axes[0].hist(es_data['cos_theta_e'], bins=50, range=(0, 1),
                alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.8)
    axes[0].set_xlabel(r'cos($\theta_{e}$) (Electron Angle)', fontsize=14)
    axes[0].set_ylabel('Counts', fontsize=14)
    axes[0].set_title('ES Events: Electron Angular Distribution', fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].axvline(0.5, color='red', linestyle='--', alpha=0.5, linewidth=2, label='60°')
    axes[0].legend()
    axes[0].text(0.05, 0.95, 'Electrons always\nforward (0° - 90°)', transform=axes[0].transAxes,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    axes[1].hist(cc_data['cos_theta'], bins=50, range=(-1, 1),
                alpha=0.7, color='coral', edgecolor='black', linewidth=0.8)
    axes[1].set_xlabel(r'cos($\theta$)', fontsize=14)
    axes[1].set_ylabel('Counts', fontsize=14)
    axes[1].set_title('CC Events: Angular Distribution', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    axes[1].axvline(0, color='red', linestyle='--', alpha=0.5, linewidth=2)
    axes[1].text(0.05, 0.95, 'More isotropic\ndistribution', transform=axes[1].transAxes,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "cos_theta_es_cc_comparison.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'cos_theta_es_cc_comparison.png'}")
    plt.close()
    
    # Combined plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(es_data['cos_theta_e'], bins=50, range=(0, 1),
           alpha=0.6, color='steelblue', edgecolor='black', linewidth=0.8,
           label=f'ES (N={len(es_data["cos_theta_e"])})', density=True)
    ax.hist(cc_data['cos_theta'], bins=50, range=(-1, 1),
           alpha=0.6, color='coral', edgecolor='black', linewidth=0.8,
           label=f'CC (N={len(cc_data["cos_theta"])})', density=True)
    
    ax.set_xlabel(r'cos($\theta_{e}$) (Electron Angle)', fontsize=14)
    ax.set_ylabel('Normalized Counts', fontsize=14)
    ax.set_title('Electron Angular Distribution: ES vs CC Events', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axvline(0, color='red', linestyle='--', alpha=0.5, linewidth=2, label='90°')
    ax.legend(fontsize=12, framealpha=0.9)
    ax.set_xlim(-1, 1)
    ax.text(0.3, 0.95, 'ES: Forward only\n(cos θ_e > 0)', transform=ax.transAxes,
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "cos_theta_combined.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'cos_theta_combined.png'}")
    plt.close()


def plot_2d_energy_angle():
    """Plot 2D energy-angle correlation with GKVM spectrum"""
    print("Generating ES events for 2D plot (GKVM)...")
    es_data = generate_es_events(50000)
    
    print("Generating CC events for 2D plot (GKVM)...")
    cc_data = generate_cc_events(50000)
    
    if cc_data is None:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # ES
    import matplotlib as mpl
    h = axes[0].hist2d(es_data['cos_theta_e'], es_data['E_e'],
                      bins=[50, 50], range=[[0, 1], [0, 70]],
                      cmap='YlOrRd', norm=mpl.colors.LogNorm())
    axes[0].set_xlabel(r'cos($\theta_{e}$) (Electron Angle)', fontsize=14)
    axes[0].set_ylabel('Electron Energy [MeV]', fontsize=14)
    axes[0].set_title('ES: Energy-Angle Correlation (GKVM Spectrum)', fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3, linestyle='--')
    cbar0 = plt.colorbar(h[3], ax=axes[0])
    cbar0.set_label('Counts', fontsize=12)
    
    # CC
    h = axes[1].hist2d(cc_data['cos_theta'], cc_data['E_e'],
                      bins=[50, 50], range=[[-1, 1], [0, 70]],
                      cmap='YlGnBu', norm=mpl.colors.LogNorm())
    axes[1].set_xlabel(r'cos($\theta_{e}$) (Electron Angle)', fontsize=14)
    axes[1].set_ylabel('Electron Energy [MeV]', fontsize=14)
    axes[1].set_title('CC: Energy-Angle Correlation (GKVM Spectrum)', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3, linestyle='--')
    cbar1 = plt.colorbar(h[3], ax=axes[1])
    cbar1.set_label('Counts', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "energy_angle_2d_correlation.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'energy_angle_2d_correlation.png'}")
    plt.close()


def plot_2d_energy_angle_flat():
    """Plot 2D energy-angle correlation with flat spectrum"""
    print("Generating ES events for 2D plot (Flat spectrum)...")
    es_data = generate_es_events_flat(50000, E_min=1, E_max=100)
    
    print("Generating CC events for 2D plot (Flat spectrum)...")
    cc_data = generate_cc_events_flat(50000, E_min=1, E_max=100)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # ES
    import matplotlib as mpl
    h = axes[0].hist2d(es_data['cos_theta_e'], es_data['E_e'],
                      bins=[50, 50], range=[[0, 1], [0, 70]],
                      cmap='YlOrRd', norm=mpl.colors.LogNorm())
    axes[0].set_xlabel(r'cos($\theta_{e}$) (Electron Angle)', fontsize=14)
    axes[0].set_ylabel('Electron Energy [MeV]', fontsize=14)
    axes[0].set_title('ES: Energy-Angle Correlation (Flat Spectrum)', fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3, linestyle='--')
    cbar0 = plt.colorbar(h[3], ax=axes[0])
    cbar0.set_label('Counts', fontsize=12)
    
    # CC
    h = axes[1].hist2d(cc_data['cos_theta'], cc_data['E_e'],
                      bins=[50, 50], range=[[-1, 1], [0, 70]],
                      cmap='YlGnBu', norm=mpl.colors.LogNorm())
    axes[1].set_xlabel(r'cos($\theta_{e}$) (Electron Angle)', fontsize=14)
    axes[1].set_ylabel('Electron Energy [MeV]', fontsize=14)
    axes[1].set_title('CC: Energy-Angle Correlation (Flat Spectrum)', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3, linestyle='--')
    cbar1 = plt.colorbar(h[3], ax=axes[1])
    cbar1.set_label('Counts', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "energy_angle_2d_correlation_flat.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'energy_angle_2d_correlation_flat.png'}")
    plt.close()


def plot_pointing_metrics():
    """Plot quantitative pointing metrics"""
    print("Generating events for pointing analysis...")
    es_data = generate_es_events(50000)
    cc_data = generate_cc_events(50000)
    
    if cc_data is None:
        return
    
    # Mean cos(theta) vs energy
    fig, ax = plt.subplots(figsize=(12, 7))
    
    energy_bins = np.linspace(5, 50, 10)
    bin_centers = (energy_bins[:-1] + energy_bins[1:]) / 2
    
    # ES
    mean_cos_es = []
    std_cos_es = []
    for i in range(len(energy_bins) - 1):
        mask = (es_data['E_e'] >= energy_bins[i]) & (es_data['E_e'] < energy_bins[i+1])
        if np.sum(mask) > 100:
            mean_cos_es.append(np.mean(es_data['cos_theta_e'][mask]))
            std_cos_es.append(np.std(es_data['cos_theta_e'][mask]) / np.sqrt(np.sum(mask)))
        else:
            mean_cos_es.append(np.nan)
            std_cos_es.append(np.nan)
    
    # CC
    mean_cos_cc = []
    std_cos_cc = []
    for i in range(len(energy_bins) - 1):
        mask = (cc_data['E_e'] >= energy_bins[i]) & (cc_data['E_e'] < energy_bins[i+1])
        if np.sum(mask) > 100:
            mean_cos_cc.append(np.mean(cc_data['cos_theta'][mask]))
            std_cos_cc.append(np.std(cc_data['cos_theta'][mask]) / np.sqrt(np.sum(mask)))
        else:
            mean_cos_cc.append(np.nan)
            std_cos_cc.append(np.nan)
    
    mean_cos_es = np.array(mean_cos_es)
    mean_cos_cc = np.array(mean_cos_cc)
    std_cos_es = np.array(std_cos_es)
    std_cos_cc = np.array(std_cos_cc)
    
    valid_es = ~np.isnan(mean_cos_es)
    valid_cc = ~np.isnan(mean_cos_cc)
    
    ax.errorbar(bin_centers[valid_es], mean_cos_es[valid_es], yerr=std_cos_es[valid_es],
               marker='o', markersize=8, linewidth=2, capsize=5,
               label='ES Events', color='steelblue', alpha=0.8)
    ax.errorbar(bin_centers[valid_cc], mean_cos_cc[valid_cc], yerr=std_cos_cc[valid_cc],
               marker='s', markersize=8, linewidth=2, capsize=5,
               label='CC Events', color='coral', alpha=0.8)
    
    ax.axhline(0, color='black', linestyle='--', alpha=0.5, linewidth=1.5, label='Isotropic')
    ax.set_xlabel('Electron Energy [MeV]', fontsize=14)
    ax.set_ylabel(r'Mean cos($\theta$)', fontsize=14)
    ax.set_title('Forward-Peaking vs Energy', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='lower right', framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.2, 1.0)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "mean_costheta_vs_energy.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'mean_costheta_vs_energy.png'}")
    plt.close()


def main():
    """Generate all physics-based plots"""
    print("=" * 70)
    print("Generating Physics-Based Plots from GKVM Flux")
    print("=" * 70)
    
    np.random.seed(42)  # For reproducibility
    
    print("\n1. Electron energy spectra...")
    plot_electron_energies()
    
    print("\n2. Angular distributions...")
    plot_cos_theta_distributions()
    
    print("\n3. 2D energy-angle correlations (GKVM)...")
    plot_2d_energy_angle()
    
    print("\n4. 2D energy-angle correlations (Flat)...")
    plot_2d_energy_angle_flat()
    
    print("\n5. Pointing metrics...")
    plot_pointing_metrics()
    
    print("\n" + "=" * 70)
    print("All plots generated successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
