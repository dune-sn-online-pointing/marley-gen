#!/usr/bin/env python3
"""
Aggregate results from parallel SNB generation jobs.
Combines all job_XXX_results.npz files into single dataset.
"""

import numpy as np
from pathlib import Path
import sys

def main():
    base_dir = Path.cwd()
    data_dir = base_dir / "data" / "snb_1000"
    
    if not data_dir.exists():
        print(f"ERROR: Data directory not found: {data_dir}")
        sys.exit(1)
    
    # Find all job result files
    job_files = sorted(data_dir.glob("job_*_results.npz"))
    
    if not job_files:
        print(f"ERROR: No job result files found in {data_dir}")
        sys.exit(1)
    
    print(f"Found {len(job_files)} job result files")
    
    all_directions = []
    all_electrons = []
    
    for job_file in job_files:
        try:
            data = np.load(job_file, allow_pickle=True)
            
            directions = data['directions']
            electrons = data['electrons']
            
            print(f"  {job_file.name}: {len(directions)} bursts")
            
            all_directions.extend(directions)
            all_electrons.extend(electrons)
            
        except Exception as e:
            print(f"  ERROR loading {job_file.name}: {e}")
    
    # Convert to arrays
    all_directions = np.array(all_directions)
    
    # Count successful extractions
    n_successful = sum(1 for e in all_electrons if e is not None)
    
    print(f"\nAggregation complete:")
    print(f"  Total bursts: {len(all_directions)}")
    print(f"  Successful electron extractions: {n_successful}")
    print(f"  Failed: {len(all_directions) - n_successful}")
    
    # Save combined results
    output_file = data_dir / "snb_1000_combined.npz"
    
    np.savez(
        output_file,
        directions=all_directions,
        electrons=np.array(all_electrons, dtype=object),
        n_bursts=len(all_directions),
        n_successful=n_successful
    )
    
    print(f"\nSaved combined results to: {output_file}")
    
    # Print statistics
    if n_successful > 0:
        n_electrons_per_burst = [len(e) for e in all_electrons if e is not None]
        print(f"\nElectrons per burst:")
        print(f"  Mean: {np.mean(n_electrons_per_burst):.1f}")
        print(f"  Min: {np.min(n_electrons_per_burst)}")
        print(f"  Max: {np.max(n_electrons_per_burst)}")

if __name__ == "__main__":
    main()
