# MARLEY Supernova Utilities

Lightweight utilities for studying MARLEY-generated (Model of Argon Reaction Low
Energy Yields) samples used in the DUNE supernova pointing work. The repository
focuses on quick inspection of elastic-scattering (ES) and charged-current (CC)
events and on producing publication-ready plots derived from those samples.

## Repository Layout

- `configs/` – MARLEY configuration files (`*.js`) for generating new samples.
- `data/` – Local MARLEY ROOT outputs (git-ignored; add your own files here).
- `interacted-fluxes/` – SNOwGLoBES flux tables used for reweighting studies.
- `notebooks/` – Scratch notebooks for exploratory work (`tester.ipynb`).
- `plots/` – Generated figures; scripts write here by default.
- `scripts/` – Python helpers:
  - `analyze_es_sample.py` – quick ES sanity plots (energy + scattering angles).
  - `create_spectrum_comparison.py` – ES/CC comparisons and spectrum studies.
  - `create_physics_plots.py`, `extract_info.py`, `extract_info_lib.py`,
    `tester.py` – legacy analysis helpers kept for reference.
- `docs/` – Empty placeholder for future project notes.

## Quick Start

```bash
# Install dependencies (adapt to your environment)
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib uproot awkward pandas scipy seaborn
```

Generate MARLEY samples with the configuration in `configs/` (example):

```bash
marley configs/flat_ES.js -n 100000 -o data/flat_ES.root
```

Produce summary plots for a sample and a neutrino-energy slice, keeping inputs
in `data/` and outputs in `plots/`:

```bash
python3 scripts/analyze_es_sample.py --input data/flat_ES.root --enu-min 40 --enu-max 45
python3 scripts/create_spectrum_comparison.py
```

## Data Notes

- Large MARLEY ROOT files are ignored via `.gitignore`; store them under
  `data/` locally.
- The scripts resolve paths relative to the repository root, so they work from
  any working directory.
- The ES utilities rely on MARLEY's "l"/"r" convention (`pdgl` = neutrino,
  `pdgr` = electron) when extracting final-state kinematics.

## Contributing

Keep new utilities under `scripts/`, document additional configuration files in
`configs/`, and prefer updating this README over adding loose status notes to
the repository.

