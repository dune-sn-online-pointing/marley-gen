#!/usr/bin/env python3
"""
Compare Electron Spectra for Different Neutrino Inputs
======================================================
Compare ES electron energies for:
1. Flat neutrino spectrum (2-70 MeV) - using MARLEY flat_ES.js
2. Mono-energetic neutrinos (45 MeV) - using MARLEY mono45 config for ES
3. GKVM supernova spectrum - using MARLEY with GKVM flux

WARNING: This script requires MARLEY-generated ROOT files:
  - flat_ES.root (from flat_ES.js config)
  - mono45_ES.root (from mono45.js modified for ES)
  - gkvm_ES.root (from GKVM flux)

If these files don't exist, you need to:
1. Run: marley flat_ES.js
2. Create mono45_ES.js (copy mono45.js, change reactions to ES.react)
3. Create gkvm_ES.js (use GKVM flux with ES.react)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "plots"
FLUX_DIR = REPO_ROOT / "interacted-fluxes"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Try to import uproot for reading ROOT files
try:
    import uproot
    HAS_UPROOT = True
except ImportError:
    HAS_UPROOT = False
    print("WARNING: uproot not available. Install with: pip install uproot")

def load_marley_es_data(root_file):
    """Load ES event data from MARLEY ROOT file
    
    Returns dict with E_nu (initial), E_e (electron KE), E_nu_final (scattered neutrino),
    and angular deviations relative to initial neutrino direction
    
    For ES events (ν + e⁻ → ν' + e⁻), MARLEY stores both final state particles.
    We calculate the scattering angle relative to the initial neutrino direction.
    """
    if not HAS_UPROOT:
        raise ImportError("uproot required to read MARLEY ROOT files")
    
    root_path = Path(root_file)
    if not root_path.is_absolute():
        root_path = (REPO_ROOT / root_path).resolve()

    if not root_path.exists():
        raise FileNotFoundError(f"MARLEY file not found: {root_path}")

    with uproot.open(root_path) as f:
        tree = f["mst"]
        
        # Get initial neutrino energy
        E_nu = tree["Ev"].array(library="np")
        
        # Prefer dedicated ES branches ("l"=neutrino, "r"=electron)
        es_branches = [
            "pdgr", "pdgl", "KEr", "KEl", "pxr", "pyr", "pzr", "pxl", "pyl", "pzl",
        ]
        if all(branch in tree.keys() for branch in es_branches):
            pdg_r = tree["pdgr"].array(library="np")
            pdg_l = tree["pdgl"].array(library="np")

            E_e = tree["KEr"].array(library="np")
            E_nu_final = tree["KEl"].array(library="np")

            pxr = tree["pxr"].array(library="np")
            pyr = tree["pyr"].array(library="np")
            pzr = tree["pzr"].array(library="np")
            pxl = tree["pxl"].array(library="np")
            pyl = tree["pyl"].array(library="np")
            pzl = tree["pzl"].array(library="np")

            p_mag_e = np.sqrt(pxr**2 + pyr**2 + pzr**2)
            p_mag_nu = np.sqrt(pxl**2 + pyl**2 + pzl**2)

            cos_theta_e = np.divide(pzr, p_mag_e, out=np.ones_like(pzr), where=p_mag_e > 0)
            cos_theta_nu = np.divide(pzl, p_mag_nu, out=np.ones_like(pzl), where=p_mag_nu > 0)

            # Basic validation: ensure we indeed have an electron and neutrino
            if not np.all(np.isin(pdg_r, [11, -11])):
                print("WARNING: Unexpected PDG codes in pdgr branch")
            if not np.all(np.isin(pdg_l, [12, -12, 14, -14, 16, -16])):
                print("WARNING: Unexpected PDG codes in pdgl branch")

        else:
            # Fallback to generic particle list (older MARLEY outputs)
            pdgp = tree["pdgp"].array(library="np")
            KEp = tree["KEp"].array(library="np")

            try:
                pxl = tree["pxl"].array(library="np")
                pyl = tree["pyl"].array(library="np")
                pzl = tree["pzl"].array(library="np")
                pxr = tree["pxr"].array(library="np")
                pyr = tree["pyr"].array(library="np")
                pzr = tree["pzr"].array(library="np")
                has_momentum = True
            except Exception:
                has_momentum = False

            E_e_list = []
            E_nu_list = []
            cos_theta_e_list = []
            cos_theta_nu_list = []

            for i in range(len(pdgp)):
                event_pdgs = pdgp[i]
                event_kes = KEp[i]

                electron_mask = event_pdgs == 11
                neutrino_mask = np.isin(event_pdgs, [12, -12, 14, -14, 16, -16])

                if np.any(electron_mask):
                    E_e_list.append(event_kes[electron_mask][0])
                    if has_momentum:
                        px, py, pz = pxr[i], pyr[i], pzr[i]
                        p_mag = np.sqrt(px**2 + py**2 + pz**2)
                        cos_theta_e_list.append(pz / p_mag if p_mag > 0 else 1.0)
                    else:
                        cos_theta_e_list.append(1.0)
                else:
                    E_e_list.append(0.0)
                    cos_theta_e_list.append(1.0)

                if np.any(neutrino_mask):
                    nu_idx = np.where(neutrino_mask)[0][0]
                    E_nu_list.append(event_kes[nu_idx])
                    if has_momentum:
                        px, py, pz = pxl[i], pyl[i], pzl[i]
                        p_mag = np.sqrt(px**2 + py**2 + pz**2)
                        cos_theta_nu_list.append(pz / p_mag if p_mag > 0 else 1.0)
                    else:
                        cos_theta_nu_list.append(1.0)
                else:
                    E_nu_list.append(E_nu[i])
                    cos_theta_nu_list.append(1.0)

            E_e = np.array(E_e_list)
            E_nu_final = np.array(E_nu_list)
            cos_theta_e = np.array(cos_theta_e_list)
            cos_theta_nu = np.array(cos_theta_nu_list)
        
    return {
        'E_nu': E_nu[:len(E_e)], 
        'E_e': E_e, 
        'E_nu_final': E_nu_final,
        'cos_theta_e': cos_theta_e,
        'cos_theta_nu': cos_theta_nu
    }

def load_marley_cc_data(root_file):
    """Load CC event data from MARLEY ROOT file
    
    For CC events (ν + nucleus → e⁻ + nucleus'), there is NO final state neutrino,
    only the produced electron.
    
    Returns dict with E_nu (initial), E_e (produced electron KE)
    """
    if not HAS_UPROOT:
        raise ImportError("uproot required to read MARLEY ROOT files")
    
    root_path = Path(root_file)
    if not root_path.is_absolute():
        root_path = (REPO_ROOT / root_path).resolve()

    if not root_path.exists():
        raise FileNotFoundError(f"MARLEY file not found: {root_path}")

    with uproot.open(root_path) as f:
        tree = f["mst"]
        
        # Get initial neutrino energy
        E_nu = tree["Ev"].array(library="np")
        
        # Final state particles
        pdgp = tree["pdgp"].array(library="np")
        KEp = tree["KEp"].array(library="np")
        
        # Extract electron from each event
        E_e = []
        for i in range(len(pdgp)):
            event_pdgs = pdgp[i]
            event_kes = KEp[i]
            
            # Find electron (PDG 11)
            electron_mask = event_pdgs == 11
            
            if np.any(electron_mask):
                E_e.append(event_kes[electron_mask][0])
            else:
                E_e.append(0.0)
        
        E_e = np.array(E_e)
        
    return {'E_nu': E_nu[:len(E_e)], 'E_e': E_e}

def es_kinematics(E_nu, cos_theta_nu):
    """Calculate electron recoil energy and final neutrino energy for ES scattering
    
    Args:
        E_nu: Initial neutrino energy (MeV)
        cos_theta_nu: Cosine of NEUTRINO scattering angle
    
    Returns:
        T_e: Electron recoil kinetic energy (MeV)
        E_nu_final: Final neutrino energy (MeV)
        
    Formulas: 
        T_e = 2·E_ν²·sin²(θ_ν/2) / (m_e + 2·E_ν·sin²(θ_ν/2))
        E_ν' = E_ν - T_e  (energy conservation)
    """
    m_e = 0.511  # MeV
    sin2_half = (1.0 - cos_theta_nu) / 2.0
    numerator = 2 * E_nu**2 * sin2_half
    denominator = m_e + 2 * E_nu * sin2_half
    T_e = numerator / denominator
    E_nu_final = E_nu - T_e  # Energy conservation
    return T_e, E_nu_final


def compute_es_electron_cos_theta(E_nu, cos_theta_nu, T_e):
    """Compute electron scattering angle cosine from kinematics.

    The initial neutrino travels along +z. We construct the post-scatter
    momenta in the scattering plane to extract the electron deflection.
    """
    E_nu_prime = E_nu - T_e
    sin_theta_nu = np.sqrt(np.clip(1.0 - cos_theta_nu**2, 0.0, 1.0))

    # Momentum components in the scattering plane
    px = -E_nu_prime * sin_theta_nu
    pz = E_nu - E_nu_prime * cos_theta_nu
    p_mag = np.sqrt(px**2 + pz**2)

    with np.errstate(divide='ignore', invalid='ignore'):
        cos_theta_e = np.where(p_mag > 0, pz / p_mag, 1.0)

    return cos_theta_e

def generate_flat_es(n_events=50000, E_min=1, E_max=100):
    """Generate ES events with flat neutrino spectrum"""
    E_nu = np.random.uniform(E_min, E_max, n_events)
    
    # Sample neutrino scattering angle from ES differential cross section
    # dσ/d(cos θ) ∝ (1 + cos θ)^2
    cos_theta_nu = []
    while len(cos_theta_nu) < n_events:
        cos_theta_proposal = np.random.uniform(-1, 1, n_events * 2)
        weight = (1 + cos_theta_proposal)**2 / 4.0
        accept = np.random.uniform(0, 1, len(cos_theta_proposal)) < weight
        cos_theta_nu.extend(cos_theta_proposal[accept])
    cos_theta_nu = np.array(cos_theta_nu[:n_events])
    
    E_e, E_nu_final = es_kinematics(E_nu, cos_theta_nu)
    cos_theta_e = compute_es_electron_cos_theta(E_nu, cos_theta_nu, E_e)
    return {
        'E_nu': E_nu,
        'E_e': E_e,
        'E_nu_final': E_nu_final,
        'cos_theta_nu': cos_theta_nu,
        'cos_theta_e': cos_theta_e,
    }

def select_quasi_mono_es(flat_data, E_nu_min=45.0, E_nu_max=50.0):
    """Select events from flat spectrum data where neutrino energy is in narrow range
    
    This gives us a quasi-monoenergetic sample from existing data
    """
    E_nu = flat_data['E_nu']
    E_e = flat_data['E_e']
    E_nu_final = flat_data['E_nu_final']
    cos_theta_e = flat_data.get('cos_theta_e', np.zeros_like(E_nu))
    cos_theta_nu = flat_data.get('cos_theta_nu', np.ones_like(E_nu))
    
    # Select events in the energy window
    mask = (E_nu >= E_nu_min) & (E_nu < E_nu_max)
    
    return {
        'E_nu': E_nu[mask],
        'E_e': E_e[mask],
        'E_nu_final': E_nu_final[mask],
        'cos_theta_e': cos_theta_e[mask],
        'cos_theta_nu': cos_theta_nu[mask],
        'E_nu_min': E_nu_min,
        'E_nu_max': E_nu_max
    }

def generate_mono_es(n_events=50000, E_nu_fixed=45.0):
    """Generate ES events with mono-energetic neutrinos (FALLBACK ONLY)"""
    E_nu = np.full(n_events, E_nu_fixed)
    
    # Sample neutrino scattering angle from ES differential cross section
    cos_theta_nu = []
    while len(cos_theta_nu) < n_events:
        cos_theta_proposal = np.random.uniform(-1, 1, n_events * 2)
        weight = (1 + cos_theta_proposal)**2 / 4.0
        accept = np.random.uniform(0, 1, len(cos_theta_proposal)) < weight
        cos_theta_nu.extend(cos_theta_proposal[accept])
    cos_theta_nu = np.array(cos_theta_nu[:n_events])
    
    E_e, E_nu_final = es_kinematics(E_nu, cos_theta_nu)
    cos_theta_e = compute_es_electron_cos_theta(E_nu, cos_theta_nu, E_e)
    return {
        'E_nu': E_nu,
        'E_e': E_e,
        'E_nu_final': E_nu_final,
        'cos_theta_nu': cos_theta_nu,
        'cos_theta_e': cos_theta_e,
    }

def generate_gkvm_es(n_events=50000):
    """Generate ES events from GKVM flux"""
    # Load GKVM flux
    flavors = ["nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
    files = [FLUX_DIR / f"gvkm_{flav}_e_ar40kt_events.dat" for flav in flavors]
    
    spectra = []
    for f in files:
        if f.exists():
            data = np.genfromtxt(f, delimiter=' ', skip_footer=2)
            spectra.append(data[:, 1])
    
    data = np.genfromtxt(files[0], delimiter=' ', skip_footer=2)
    E_nu_mev = data[:, 0] * 1000
    flux = np.sum(spectra, axis=0)
    flux_norm = flux / np.sum(flux)
    
    E_nu = np.random.choice(E_nu_mev, size=n_events, p=flux_norm)
    
    # Sample neutrino scattering angle from ES differential cross section
    cos_theta_nu = []
    while len(cos_theta_nu) < n_events:
        cos_theta_proposal = np.random.uniform(-1, 1, n_events * 2)
        weight = (1 + cos_theta_proposal)**2 / 4.0
        accept = np.random.uniform(0, 1, len(cos_theta_proposal)) < weight
        cos_theta_nu.extend(cos_theta_proposal[accept])
    cos_theta_nu = np.array(cos_theta_nu[:n_events])
    
    E_e, E_nu_final = es_kinematics(E_nu, cos_theta_nu)
    cos_theta_e = compute_es_electron_cos_theta(E_nu, cos_theta_nu, E_e)
    return {
        'E_nu': E_nu,
        'E_e': E_e,
        'E_nu_final': E_nu_final,
        'cos_theta_nu': cos_theta_nu,
        'cos_theta_e': cos_theta_e,
    }

def plot_es_comparison():
    """Create 3-panel comparison of ES electron energies"""
    
    # Try to load MARLEY data first, fall back to manual generation
    print("\n" + "="*70)
    print("IMPORTANT: This script should use MARLEY-generated ES samples!")
    print("="*70)
    
    # Check for MARLEY ROOT files
    marley_files = {
        'flat': 'data/flat_ES.root',
        'mono': 'data/mono45_ES.root',
        'gkvm': 'data/gkvm_ES.root'
    }
    
    use_marley = all(Path(f).exists() for f in marley_files.values()) and HAS_UPROOT
    
    if use_marley:
        print("✓ Loading MARLEY-generated ES samples...")
        flat_data = load_marley_es_data(marley_files['flat'])
        gkvm_data = load_marley_es_data(marley_files['gkvm'])
    else:
        print("✗ MARLEY files not found. Using manual ES kinematics (NOT RECOMMENDED)")
        print("  Generate proper samples with MARLEY for accurate results!")
        print("  Missing files:")
        for name, path in marley_files.items():
            if not Path(path).exists():
                print(f"    - {path}")
        print()
        
        print("Generating flat spectrum ES (manual)...")
        flat_data = generate_flat_es(50000)
        
        print("Generating GKVM spectrum ES (manual)...")
        gkvm_data = generate_gkvm_es(50000)
    
    # Select quasi-monoenergetic sample from flat spectrum data
    print("Selecting quasi-mono sample (45-50 MeV) from flat spectrum...")
    mono_data = select_quasi_mono_es(flat_data, E_nu_min=45.0, E_nu_max=50.0)
    print(f"  Selected {len(mono_data['E_nu'])} events")
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Flat spectrum
    axes[0].hist(flat_data['E_e'], bins=100, range=(0, 60),
                alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
    axes[0].set_xlabel('Electron Energy [MeV]', fontsize=13)
    axes[0].set_ylabel('Counts', fontsize=13)
    axes[0].set_title('ES: Flat Neutrino Spectrum\n(ν energy: 1-100 MeV)', 
                     fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log')
    axes[0].text(0.05, 0.95, f'Mean E_e = {np.mean(flat_data["E_e"]):.2f} MeV',
                transform=axes[0].transAxes, fontsize=11,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    # Quasi-monoenergetic (selected from flat spectrum)
    # Show initial neutrino, final neutrino (ES only), and final electron energies
    axes[1].hist(mono_data['E_nu'], bins=50, range=(0, 60),
                alpha=0.5, color='blue', edgecolor='black', linewidth=0.8,
                label=f'Initial ν ({mono_data["E_nu_min"]:.0f}-{mono_data["E_nu_max"]:.0f} MeV)')
    axes[1].hist(mono_data['E_nu_final'], bins=50, range=(0, 60),
                alpha=0.5, color='red', edgecolor='black', linewidth=0.8,
                label='Final ν (scattered)')
    axes[1].hist(mono_data['E_e'], bins=50, range=(0, 60),
                alpha=0.5, color='orange', edgecolor='black', linewidth=0.8,
                label='Final e⁻ (recoil)')
    axes[1].set_xlabel('Energy [MeV]', fontsize=13)
    axes[1].set_ylabel('Counts', fontsize=13)
    axes[1].set_title(f'ES: Quasi-Monoenergetic Neutrinos\n(Initial ν: {mono_data["E_nu_min"]:.0f}-{mono_data["E_nu_max"]:.0f} MeV)', 
                     fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=9, loc='upper right')
    axes[1].text(0.05, 0.95, 
                f'N events = {len(mono_data["E_nu"])}\n' +
                f'Initial: ν={np.mean(mono_data["E_nu"]):.2f} MeV\n' +
                f'Final: ν={np.mean(mono_data["E_nu_final"]):.2f} MeV, e⁻={np.mean(mono_data["E_e"]):.2f} MeV',
                transform=axes[1].transAxes, fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # GKVM spectrum
    axes[2].hist(gkvm_data['E_e'], bins=100, range=(0, 60),
                alpha=0.7, color='forestgreen', edgecolor='black', linewidth=0.5)
    axes[2].set_xlabel('Electron Energy [MeV]', fontsize=13)
    axes[2].set_ylabel('Counts', fontsize=13)
    axes[2].set_title('ES: GKVM SN Spectrum\n(Garching Model)', 
                     fontsize=13, fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    axes[2].set_yscale('log')
    axes[2].text(0.05, 0.95, f'Mean E_e = {np.mean(gkvm_data["E_e"]):.2f} MeV',
                transform=axes[2].transAxes, fontsize=11,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "es_electron_energy_comparison.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'es_electron_energy_comparison.png'}")
    plt.close()
    
    # Also create overlay comparison
    fig, ax = plt.subplots(figsize=(12, 7))
    
    ax.hist(flat_data['E_e'], bins=80, range=(0, 60),
           alpha=0.5, color='steelblue', edgecolor='black', linewidth=0.8,
           label='Flat ν spectrum (1-100 MeV)', density=True)
    ax.hist(mono_data['E_e'], bins=80, range=(0, 60),
           alpha=0.5, color='orange', edgecolor='black', linewidth=0.8,
           label=f'Quasi-mono ν ({mono_data["E_nu_min"]:.0f}-{mono_data["E_nu_max"]:.0f} MeV)', density=True)
    ax.hist(gkvm_data['E_e'], bins=80, range=(0, 60),
           alpha=0.5, color='forestgreen', edgecolor='black', linewidth=0.8,
           label='GKVM SN spectrum', density=True)
    
    ax.set_xlabel('Electron Energy [MeV]', fontsize=14)
    ax.set_ylabel('Normalized Counts', fontsize=14)
    ax.set_title('ES Electron Energy: Comparison of Neutrino Spectra', 
                fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, framealpha=0.9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "es_electron_energy_overlay.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'es_electron_energy_overlay.png'}")
    plt.close()

def plot_neutrino_spectra():
    """Show the input neutrino spectra being compared"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Flat
    E_flat = np.linspace(1, 100, 100)
    flux_flat = np.ones_like(E_flat)
    ax.plot(E_flat, flux_flat / np.max(flux_flat), 
           linewidth=3, label='Flat spectrum', color='steelblue')
    
    # Mono-energetic (delta function represented as narrow Gaussian)
    E_mono = np.linspace(1, 100, 1000)
    flux_mono = np.exp(-((E_mono - 45)**2) / (2 * 0.5**2))
    ax.plot(E_mono, flux_mono / np.max(flux_mono), 
           linewidth=3, label='Mono-energetic (45 MeV)', color='orange')
    
    # GKVM
    flavors = ["nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
    files = [FLUX_DIR / f"gvkm_{flav}_e_ar40kt_events.dat" for flav in flavors]
    
    spectra = []
    for f in files:
        if f.exists():
            data = np.genfromtxt(f, delimiter=' ', skip_footer=2)
            spectra.append(data[:, 1])
    
    data = np.genfromtxt(files[0], delimiter=' ', skip_footer=2)
    E_gkvm = data[:, 0] * 1000  # Convert to MeV
    flux_gkvm = np.sum(spectra, axis=0)
    ax.plot(E_gkvm, flux_gkvm / np.max(flux_gkvm), 
           linewidth=3, label='GKVM SN spectrum', color='forestgreen')
    
    ax.set_xlabel('Neutrino Energy [MeV]', fontsize=14)
    ax.set_ylabel('Relative Flux (normalized)', fontsize=14)
    ax.set_title('Input Neutrino Energy Spectra', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.1)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "input_neutrino_spectra.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'input_neutrino_spectra.png'}")
    plt.close()

def plot_es_vs_cc_comparison():
    """Compare ES and CC electron energy spectra for quasi-mono neutrinos"""
    print("\n3. ES vs CC electron energy comparison...")
    
    # Check for MARLEY files
    es_file = 'data/flat_ES.root'
    cc_file = 'data/flat.root'
    
    use_marley = Path(es_file).exists() and Path(cc_file).exists() and HAS_UPROOT
    
    if use_marley:
        print("✓ Loading MARLEY ES and CC samples...")
        es_flat = load_marley_es_data(es_file)
        cc_flat = load_marley_cc_data(cc_file)
    else:
        print("✗ MARLEY files not found. Using manual generation...")
        es_flat = generate_flat_es(50000)
        # For CC, we'll just use the neutrino energy as electron energy (approximation)
        # This is wrong physics but just for demonstration
        print("  WARNING: CC generation not properly implemented in fallback!")
        cc_nu = np.random.uniform(1, 100, 50000)
        cc_flat = {'E_nu': cc_nu, 'E_e': cc_nu * 0.8}  # Rough approximation
    
    # Select quasi-mono samples
    es_mono = select_quasi_mono_es(es_flat, E_nu_min=45.0, E_nu_max=50.0)
    
    # For CC, select same neutrino energy range
    mask_cc = (cc_flat['E_nu'] >= 45.0) & (cc_flat['E_nu'] < 50.0)
    cc_mono = {
        'E_nu': cc_flat['E_nu'][mask_cc],
        'E_e': cc_flat['E_e'][mask_cc]
    }
    
    print(f"  ES: {len(es_mono['E_nu'])} events")
    print(f"  CC: {len(cc_mono['E_nu'])} events")
    
    # Create comparison plot with 4 panels (2x2 grid)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # ES energy panel (top left)
    axes[0, 0].hist(es_mono['E_nu'], bins=50, range=(0, 60),
                alpha=0.5, color='blue', edgecolor='black', linewidth=0.8,
                label='Initial ν')
    axes[0, 0].hist(es_mono['E_nu_final'], bins=50, range=(0, 60),
                alpha=0.5, color='red', edgecolor='black', linewidth=0.8,
                label='Final ν')
    axes[0, 0].hist(es_mono['E_e'], bins=50, range=(0, 60),
                alpha=0.5, color='orange', edgecolor='black', linewidth=0.8,
                label='Final e⁻')
    axes[0, 0].set_xlabel('Energy [MeV]', fontsize=13)
    axes[0, 0].set_ylabel('Counts', fontsize=13)
    axes[0, 0].set_title('ES Energy\nν + e⁻ → ν + e⁻', 
                     fontsize=13, fontweight='bold')
    axes[0, 0].legend(fontsize=10)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].text(0.05, 0.95, 
                f'ν: 45-50 MeV\n' +
                f'Mean e⁻: {np.mean(es_mono["E_e"]):.1f} MeV',
                transform=axes[0, 0].transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    # ES electron angle (top right)
    axes[0, 1].hist(es_mono['cos_theta_e'], bins=50, range=(-1, 1),
                alpha=0.7, color='orange', edgecolor='black', linewidth=0.8)
    axes[0, 1].set_xlabel('cos(θ_scatter)', fontsize=13)
    axes[0, 1].set_ylabel('Counts', fontsize=13)
    axes[0, 1].set_title('ES: Electron Scattering Angle\n(relative to initial ν)', 
                     fontsize=13, fontweight='bold')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].axvline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='90°')
    axes[0, 1].axvline(1, color='blue', linestyle='--', linewidth=1, alpha=0.5, label='0° (forward)')
    axes[0, 1].legend(fontsize=9)
    axes[0, 1].text(0.05, 0.95, 
                f'Mean: {np.mean(es_mono["cos_theta_e"]):.2f}\n' +
                f'Forward-peaked',
                transform=axes[0, 1].transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # ES neutrino angle (bottom left)
    axes[1, 0].hist(es_mono['cos_theta_nu'], bins=50, range=(-1, 1),
                alpha=0.7, color='red', edgecolor='black', linewidth=0.8)
    axes[1, 0].set_xlabel('cos(θ_scatter)', fontsize=13)
    axes[1, 0].set_ylabel('Counts', fontsize=13)
    axes[1, 0].set_title('ES: Neutrino Scattering Angle\n(relative to initial ν)', 
                     fontsize=13, fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axvline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='90°')
    axes[1, 0].axvline(1, color='blue', linestyle='--', linewidth=1, alpha=0.5, label='0° (forward)')
    axes[1, 0].legend(fontsize=9)
    axes[1, 0].text(0.05, 0.95, 
                f'Mean: {np.mean(es_mono["cos_theta_nu"]):.2f}\n' +
                f'Strongly forward',
                transform=axes[1, 0].transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
    
    # CC energy panel (bottom right)
    axes[1, 1].hist(cc_mono['E_nu'], bins=50, range=(0, 60),
                alpha=0.5, color='blue', edgecolor='black', linewidth=0.8,
                label='Initial ν')
    axes[1, 1].hist(cc_mono['E_e'], bins=50, range=(0, 60),
                alpha=0.5, color='orange', edgecolor='black', linewidth=0.8,
                label='Final e⁻')
    axes[1, 1].set_xlabel('Energy [MeV]', fontsize=13)
    axes[1, 1].set_ylabel('Counts', fontsize=13)
    axes[1, 1].set_title('CC Energy\nν + ⁴⁰Ar → e⁻ + ⁴⁰K*', 
                     fontsize=13, fontweight='bold')
    axes[1, 1].legend(fontsize=10)
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].text(0.05, 0.95, 
                f'ν: 45-50 MeV\n' +
                f'Mean e⁻: {np.mean(cc_mono["E_e"]):.1f} MeV\n' +
                f'NO final ν!',
                transform=axes[1, 1].transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "es_vs_cc_comparison.png", dpi=300)
    print(f"Saved: {OUT_DIR / 'es_vs_cc_comparison.png'}")
    plt.close()

def main():
    """Generate comparison plots"""
    print("=" * 70)
    print("Creating ES Electron Energy Comparisons")
    print("=" * 70)
    
    np.random.seed(42)
    
    print("\n1. Input neutrino spectra...")
    plot_neutrino_spectra()
    
    print("\n2. ES electron energy comparison...")
    plot_es_comparison()
    
    print("\n3. ES vs CC comparison (quasi-mono neutrinos)...")
    plot_es_vs_cc_comparison()
    
    print("\n" + "=" * 70)
    print("Comparison plots complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
