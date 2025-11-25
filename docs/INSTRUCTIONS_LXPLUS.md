# Running the Toy Study on LXPLUS

## Overview
This toy study generates 1000 MARLEY samples (each with 400 ES electrons from a random neutrino direction) and reconstructs each direction using MCMC (emcee with 64 walkers, 2000 steps).

**Estimated Runtime:** ~17-28 hours for 1000 samples (10-15 seconds per sample)

## Prerequisites on LXPLUS

### 1. Set up MARLEY
```bash
# If MARLEY is not installed, download and build it
cd ~/
git clone https://github.com/MARLEY-MC/marley.git marley-1.2.0
cd marley-1.2.0
mkdir build && cd build

# Load dependencies
module load gsl/2.7
module load root/6.28.04

# Build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/marley-1.2.0/install
make -j4
make install

# Set environment
export MARLEY=$HOME/marley-1.2.0
export PATH=$MARLEY/build:$PATH
export LD_LIBRARY_PATH=$MARLEY/build:$LD_LIBRARY_PATH
```

### 2. Set up Python environment
```bash
# On lxplus, use Python 3.9+
module load python/3.9

# Install required packages
pip install --user numpy matplotlib emcee tqdm scipy uproot awkward
```

## Running the Toy Study

### Option 1: Direct execution (for testing small samples)
```bash
cd marley-gen
python3 scripts/toy_study_with_checkpoints.py
```

### Option 2: Submit to HTCondor (RECOMMENDED for 1000 samples)

Create a HTCondor submission script `submit_toy_study.sub`:
```
# submit_toy_study.sub
executable = /usr/bin/python3
arguments = scripts/toy_study_with_checkpoints.py
output = logs/toy_study_$(ClusterId).out
error = logs/toy_study_$(ClusterId).err
log = logs/toy_study_$(ClusterId).log

# Request resources
request_cpus = 1
request_memory = 4GB
request_disk = 2GB

# Environment
+MaxRuntime = 172800  # 48 hours
environment = "MARLEY=/afs/cern.ch/user/YOUR_USERNAME/marley-1.2.0 LD_LIBRARY_PATH=/afs/cern.ch/user/YOUR_USERNAME/marley-1.2.0/build:$LD_LIBRARY_PATH PATH=/afs/cern.ch/user/YOUR_USERNAME/marley-1.2.0/build:$PATH"

# Transfer files
should_transfer_files = YES
transfer_input_files = configs/, scripts/toy_study_with_checkpoints.py
transfer_output_files = data/toy_study/reconstruction_results.npz, data/toy_study/true_directions.npy, plots/toy_study_1000samples_mcmc_reconstruction.png
when_to_transfer_output = ON_EXIT

queue 1
```

Then submit:
```bash
mkdir -p logs
condor_submit submit_toy_study.sub
condor_q  # Check status
```

### Option 3: Parallel processing with job array

For faster execution, split into multiple jobs. Modify the script to accept start/end indices:

Create `scripts/toy_study_parallel.py` (see repository) and submit array job:
```
# submit_parallel.sub
executable = /usr/bin/python3
arguments = scripts/toy_study_parallel.py $(Process) 100  # Process 100 samples per job
output = logs/toy_$(Process).out
error = logs/toy_$(Process).err
log = logs/toy_$(Process).log

request_cpus = 1
request_memory = 4GB
request_disk = 1GB
+MaxRuntime = 86400  # 24 hours

queue 10  # 10 jobs × 100 samples = 1000 samples
```

Then aggregate results:
```bash
python3 scripts/aggregate_parallel_results.py
```

## Files Generated

- `data/toy_study/true_directions.npy` - True neutrino directions for each sample
- `data/toy_study/reconstruction_results.npz` - MCMC reconstruction results
- `data/toy_study/checkpoint.npz` - Checkpoint file (for resuming interrupted runs)
- `plots/toy_study_1000samples_mcmc_reconstruction.png` - Final plot with 68% quantile

## Expected Output

The final plot will show:
- Distribution of cos(θ_true - θ_reconstructed) for 1000 random directions
- 68% quantile line (red dashed) showing reconstruction precision
- Mean reconstruction accuracy (green dashed)
- Statistics box with performance metrics

Typical results should show:
- Mean cos(θ): ~0.95-0.98 (very good reconstruction)
- 68% precision: ~10-20 degrees
- Most reconstructions within a few degrees of true direction

## Troubleshooting

### MARLEY not found
```bash
export MARLEY=$HOME/marley-1.2.0
export PATH=$MARLEY/build:$PATH
export LD_LIBRARY_PATH=$MARLEY/build:$LD_LIBRARY_PATH
```

### Python packages missing
```bash
pip install --user numpy matplotlib emcee tqdm scipy uproot awkward
```

### ROOT file extraction fails
Check that MARLEY was built with ROOT support:
```bash
marley --version  # Should show "ROOT support: yes"
```

### Jobs running too long
Consider splitting into smaller batches (Option 3) or reducing:
- `n_events=400` → `n_events=200` (fewer electrons per sample)
- `nsteps=2000` → `nsteps=1000` (faster MCMC)
- `nwalkers=64` → `nwalkers=32` (fewer walkers)

## Contact
For questions about this study, see the main README.md or check the repository history.
