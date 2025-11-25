# Generating 1000 SNB Events with MARLEY on HTCondor

## Overview
Generate 1000 supernova burst (SNB) events, each with:
- 400 ES (elastic scattering) interactions
- Same random neutrino direction per burst
- Different directions for each burst

Split across **100 HTCondor jobs** (10 bursts per job) for parallel execution.

## Quick Start

### 1. Prepare directories
```bash
cd /afs/cern.ch/work/e/evilla/private/dune/marley-gen
mkdir -p logs data/snb_1000
```

### 2. Submit jobs to HTCondor
```bash
condor_submit condor_snb_1000.sub
```

### 3. Monitor jobs
```bash
# Check job status
condor_q

# Watch progress (optional)
watch -n 10 'condor_q; ls -1 data/snb_1000/job_*_results.npz | wc -l'
```

### 4. Wait for completion
All 100 jobs should complete in **30-60 minutes** (depending on cluster load).

### 5. Aggregate results
```bash
python3 scripts/aggregate_snb_results.py
```

This creates: `data/snb_1000/snb_1000_combined.npz`

## Output Files

### Individual job results
- `data/snb_1000/job_000_results.npz` through `job_099_results.npz`
- Each contains:
  - `directions`: (10, 3) array of neutrino directions
  - `electrons`: list of electron momentum arrays (N×3)
  - `job_id`, `n_bursts`: metadata

### Combined results
- `data/snb_1000/snb_1000_combined.npz`
  - `directions`: (1000, 3) array of all neutrino directions
  - `electrons`: list of 1000 electron arrays
  - `n_bursts`: total count
  - `n_successful`: successful electron extractions

## File Structure
```
marley-gen/
├── condor_snb_1000.sub          # HTCondor submission script
├── scripts/
│   ├── generate_snb_batch.py    # Per-job SNB generator
│   └── aggregate_snb_results.py # Combine all results
├── data/snb_1000/
│   ├── job_000_results.npz      # Job 0 output (10 bursts)
│   ├── job_001_results.npz      # Job 1 output (10 bursts)
│   ├── ...
│   ├── job_099_results.npz      # Job 99 output (10 bursts)
│   └── snb_1000_combined.npz    # Final combined dataset
└── logs/
    ├── snb_0.out, snb_0.err     # Job 0 logs
    ├── snb_1.out, snb_1.err     # Job 1 logs
    └── ...
```

## Checking Progress

### Number of completed jobs
```bash
ls -1 data/snb_1000/job_*_results.npz | wc -l
```

### Check specific job output
```bash
cat logs/snb_0.out   # Check job 0 stdout
cat logs/snb_0.err   # Check job 0 stderr
```

### Check for errors
```bash
grep -l ERROR logs/snb_*.err
```

## Troubleshooting

### Jobs held/failed
```bash
condor_q -analyze <job_id>
condor_q -better-analyze <job_id>
```

### Resubmit specific jobs
If job 42 failed:
```bash
# Edit condor_snb_1000.sub and change last line to:
# queue 1 Process in (42)

condor_submit condor_snb_1000.sub
```

### Missing MARLEY
Ensure MARLEY is installed:
```bash
export MARLEY=$HOME/marley-1.2.0
export PATH=$MARLEY/build:$PATH
export LD_LIBRARY_PATH=$MARLEY/build:$LD_LIBRARY_PATH
marley --version
```

### Python packages missing
```bash
pip install --user numpy uproot awkward
```

## Next Steps

After generating the 1000 SNB events, you can:

1. **Run MCMC reconstruction** on each burst to recover neutrino directions
2. **Plot angular resolution** (cos θ distribution)
3. **Analyze systematic biases** in directional reconstruction
4. **Compare different reconstruction methods**

See `scripts/reconstruct_snb_directions.py` (to be created) for reconstruction pipeline.

## Performance Notes

- **Time per burst**: ~2-3 minutes (400 events + electron extraction)
- **Time per job** (10 bursts): ~20-30 minutes
- **Total walltime**: 30-60 minutes (with 100 parallel jobs)
- **Disk space**: ~50-100 MB per job, ~5-10 GB total
- **Memory**: ~1-2 GB per job

## Dataset Description

Each SNB event represents:
- A supernova burst from a **random direction** on the sky
- 400 electron-scattering interactions in liquid argon
- Electron momenta extracted for directional reconstruction

This dataset can be used to study:
- Supernova neutrino directionality in DUNE
- Angular resolution vs. statistics (400 events per burst)
- Systematic effects in directional reconstruction
- Background rejection strategies
