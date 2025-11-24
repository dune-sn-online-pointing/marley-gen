# Quick Start on LXPLUS

## What's Ready
✅ Generated 400 ES electrons from direction (1,1,1) - **DONE**
✅ Plotted electron angular distributions - **DONE** 
✅ Validated cos(θ_ν - θ_e) with 68% quantile - **DONE**
✅ Created MCMC reconstruction script with emcee - **READY**

## What to Run on LXPLUS

### Step 1: Clone and setup
```bash
ssh lxplus.cern.ch
git clone https://github.com/dune-sn-online-pointing/marley-gen.git
cd marley-gen
```

### Step 2: Run the toy study
```bash
# Check INSTRUCTIONS_LXPLUS.md for full details
python3 scripts/toy_study_1000samples.py
```

This will:
- Generate 1000 MARLEY samples (random neutrino directions)
- Each sample: 400 ES electrons  
- Run MCMC reconstruction (64 walkers, 2000 steps) per sample
- Output: `plots/toy_study_1000samples_mcmc_reconstruction.png`

**Time:** ~17-28 hours (can be split into parallel jobs on HTCondor)

## Key Files
- `scripts/toy_study_1000samples.py` - Main script
- `INSTRUCTIONS_LXPLUS.md` - Detailed instructions
- `configs/directional_ES.js` - Example MARLEY config

## Already Generated (for reference)
- `plots/directional_400_comprehensive.png` - Electron angles from direction (1,1,1)
- `plots/cos_angle_separation_68percentile.png` - Angular separation with 68% quantile
- `data/directional_electron_data.txt` - Extracted electron momenta

## Expected Output
Final plot showing:
- cos(θ_true - θ_reconstructed) distribution for 1000 samples
- Red line: 68% quantile (reconstruction precision)
- Green line: Mean reconstruction accuracy
- Statistics: typically ~15-20° precision at 68% level
