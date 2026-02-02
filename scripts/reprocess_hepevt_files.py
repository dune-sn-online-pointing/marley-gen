#!/usr/bin/python3
"""
Reprocess existing HEPEVT files to extract electron data correctly.
"""

import numpy as np
from pathlib import Path
import sys

def extract_electrons_from_hepevt(hepevt_file):
    """Extract electron momenta from MARLEY HEPEVT output."""
    try:
        electrons = []
        
        with open(hepevt_file, 'r') as f:
            lines = f.readlines()
            
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1
                continue
            
            # First line of event: event_number n_particles
            parts = line.split()
            if len(parts) >= 2:
                try:
                    n_particles = int(parts[1])
                except ValueError:
                    i += 1
                    continue
                
                # Read particle lines
                for j in range(n_particles):
                    i += 1
                    if i >= len(lines):
                        break
                    pline = lines[i].strip().split()
                    if len(pline) >= 10:
                        status = int(pline[0])
                        pdg = int(pline[1])
                        # Final state electron (status = 1, PDG = 11)
                        if status == 1 and pdg == 11:
                            px = float(pline[6])
                            py = float(pline[7])
                            pz = float(pline[8])
                            electrons.append([px, py, pz])
            i += 1
        
        return np.array(electrons) if electrons else None
    except Exception as e:
        print(f"Error extracting electrons: {e}")
        return None

def main():
    data_dir = Path("/afs/cern.ch/work/e/evilla/private/dune/marley-gen/data/snb_1000")
    
    # Process each job (0-99), 10 bursts per job
    for job_id in range(100):
        all_directions = []
        all_electrons = []
        
        for burst_idx in range(10):
            global_idx = job_id * 10 + burst_idx
            hepevt_file = data_dir / f"snb_{global_idx:04d}.hepevt"
            
            if not hepevt_file.exists():
                print(f"Missing: {hepevt_file}")
                continue
            
            # Extract direction from old results if available
            old_results = data_dir / f"job_{job_id:03d}_results.npz"
            if old_results.exists():
                old_data = np.load(old_results, allow_pickle=True)
                direction = old_data['directions'][burst_idx]
            else:
                # Generate random direction (this shouldn't happen)
                phi = np.random.uniform(0, 2*np.pi)
                cos_theta = np.random.uniform(-1, 1)
                sin_theta = np.sqrt(1 - cos_theta**2)
                x = sin_theta * np.cos(phi)
                y = sin_theta * np.sin(phi)
                z = cos_theta
                direction = np.array([x, y, z])
            
            all_directions.append(direction)
            
            # Extract electrons
            electrons = extract_electrons_from_hepevt(hepevt_file)
            all_electrons.append(electrons)
        
        # Save results
        output_file = data_dir / f"job_{job_id:03d}_results.npz"
        np.savez(
            output_file,
            directions=np.array(all_directions),
            electrons=np.array(all_electrons, dtype=object),
            job_id=job_id,
            n_bursts=10
        )
        
        successful = sum(1 for e in all_electrons if e is not None and len(e) > 0)
        print(f"Job {job_id}: {successful}/10 bursts successful")

if __name__ == "__main__":
    main()
