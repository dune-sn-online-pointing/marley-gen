#!/usr/bin/env python3
"""Quick analysis of a single MARLEY ES ROOT sample.

Given an input MARLEY elastic-scattering file this script extracts the
final-state electron ("r") and neutrino ("l") kinematics and produces
summary plots for:
  * the scattering angle cosines relative to the incoming neutrino
  * the kinetic-energy distributions of the final-state particles

Usage:
    python3 scripts/analyze_es_sample.py --input data/flat_ES.root
"""
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt

from create_spectrum_comparison import load_marley_es_data, OUT_DIR, REPO_ROOT


def make_angle_plot(es_data, output_dir: Path, tag: str) -> Path:
    """Create a histogram of scattering angle cosines for e and nu."""
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(
        es_data["cos_theta_e"], bins=60, range=(-1, 1),
        alpha=0.6, color="orange", edgecolor="black", label="Electron (r)",
    )
    ax.hist(
        es_data["cos_theta_nu"], bins=60, range=(-1, 1),
        alpha=0.6, color="royalblue", edgecolor="black", label="Neutrino (l)",
    )

    ax.set_xlabel("cos(θ_scatter)")
    ax.set_ylabel("Counts")
    ax.set_title("ES Scattering Angle Cosines (relative to initial ν)")
    ax.axvline(1, linestyle="--", color="gray", linewidth=1, alpha=0.7)
    ax.axvline(0, linestyle="--", color="gray", linewidth=1, alpha=0.7)
    ax.legend()
    ax.grid(alpha=0.3)

    text = (
        f"⟨cosθ_e⟩ = {np.mean(es_data['cos_theta_e']):.3f}\n"
        f"⟨cosθ_ν⟩ = {np.mean(es_data['cos_theta_nu']):.3f}"
    )
    ax.text(
        0.02, 0.98, text, transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{tag}_cos_theta.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def make_energy_plot(es_data, output_dir: Path, tag: str, window_label: str) -> Path:
    """Create energy spectra histograms for the final-state particles."""
    if len(es_data["E_e"]) == 0:
        raise ValueError("No ES events available after selection")

    fig, ax = plt.subplots(figsize=(8, 4))

    ax.hist(
        es_data["E_e"], bins=50, range=(0, es_data["E_e"].max() * 1.05),
        alpha=0.6, color="orange", edgecolor="black", label="Electron KE (r)",
    )
    ax.hist(
        es_data["E_nu_final"], bins=50, range=(0, es_data["E_nu_final"].max() * 1.05),
        alpha=0.6, color="royalblue", edgecolor="black", label="Neutrino KE (l)",
    )

    ax.set_xlabel("Energy [MeV]")
    ax.set_ylabel("Counts")
    # y range expand to max*1.2
    ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
    ax.set_title(f"ES Final-State Energies ({window_label})")
    ax.legend()
    ax.grid(alpha=0.3)

    text = (
        f"⟨E_e⟩ = {np.mean(es_data['E_e']):.2f} MeV\n"
        f"⟨E_ν'⟩ = {np.mean(es_data['E_nu_final']):.2f} MeV"
    )
    ax.text(
        0.02, 0.98, text, transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{tag}_energies.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze a MARLEY ES sample")
    parser.add_argument(
        "--input", "-i", type=Path,
        default=Path("data/flat_ES.root"),
        help="Path to MARLEY ES ROOT file",
    )
    parser.add_argument(
        "--tag", "-t", default="flat_ES", help="Tag used in output filenames",
    )
    parser.add_argument(
        "--enu-min", type=float, default=None,
        help="Lower bound on initial neutrino energy (MeV) for event selection",
    )
    parser.add_argument(
        "--enu-max", type=float, default=None,
        help="Upper bound on initial neutrino energy (MeV) for event selection",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.is_absolute():
        input_path = (REPO_ROOT / input_path).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"Loading ES data from {input_path} ...")
    es_data = load_marley_es_data(input_path)

    # Optional selection on initial neutrino energy
    if args.enu_min is not None or args.enu_max is not None:
        enu = es_data['E_nu']
        mask = np.ones_like(enu, dtype=bool)
        if args.enu_min is not None:
            mask &= enu >= args.enu_min
        if args.enu_max is not None:
            mask &= enu < args.enu_max

        before = len(enu)
        es_data = {key: val[mask] for key, val in es_data.items()}
        after = len(es_data['E_nu'])
        print(
            f"Applied E_nu selection: "
            f"[{args.enu_min if args.enu_min is not None else '-inf'}, "
            f"{args.enu_max if args.enu_max is not None else '+inf'}) MeV"
        )
        print(f"  Kept {after} / {before} events")

    if args.enu_min is None and args.enu_max is None:
        window_label = "all initial Eν"
    else:
        low = args.enu_min if args.enu_min is not None else "-inf"
        high = args.enu_max if args.enu_max is not None else "+inf"
        window_label = f"Eν ∈ [{low}, {high}) MeV"

    if len(es_data['E_nu']) == 0:
        raise ValueError("No ES events left after applying the selection")

    print("Creating plots...")
    angle_plot = make_angle_plot(es_data, OUT_DIR, args.tag)
    energy_plot = make_energy_plot(es_data, OUT_DIR, args.tag, window_label)

    print("Summary statistics:")
    print(f"  Events: {len(es_data['E_nu'])}")
    print(f"  Mean electron KE: {np.mean(es_data['E_e']):.2f} MeV")
    print(f"  Mean scattered ν KE: {np.mean(es_data['E_nu_final']):.2f} MeV")
    print(f"  Mean cosθ_e: {np.mean(es_data['cos_theta_e']):.3f}")
    print(f"  Mean cosθ_ν: {np.mean(es_data['cos_theta_nu']):.3f}")

    print(f"Saved plots:\n  - {angle_plot}\n  - {energy_plot}")


if __name__ == "__main__":
    main()
