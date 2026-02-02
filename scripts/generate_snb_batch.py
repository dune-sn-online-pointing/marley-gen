#!/usr/bin/python3
"""
Generate multiple SNB events with MARLEY
Each event: same random neutrino direction, 400 ES interactions

Usage: python generate_snb_batch.py <job_id> <n_bursts>
Example: python generate_snb_batch.py 0 10
"""

import sys
import os
import subprocess
import numpy as np
from pathlib import Path
import json

def generate_random_direction():
    """Generate a random unit vector (isotropic on sphere)."""
    phi = np.random.uniform(0, 2*np.pi)
    cos_theta = np.random.uniform(-1, 1)
    sin_theta = np.sqrt(1 - cos_theta**2)
    
    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta
    
    return np.array([x, y, z])

def create_marley_config(direction, output_file, n_events=400):
    """Create MARLEY JS config for a specific direction."""
    seed = int(np.random.randint(1, 2**31))
    
    config_text = f"""{{
    seed: {seed},
    
    direction: {{ x: {direction[0]}, y: {direction[1]}, z: {direction[2]} }},
    
    target: {{
      nuclides: [ 1000180400 ],
      atom_fractions: [ 1.0 ],
    }},
    
    reactions: [ "ES.react" ],
    
    source: {{
       type: "fermi-dirac",
       neutrino: "ve",
       Emin: 5,
       Emax: 50,
       temperature: 3.5,
       eta: 0
     }},
    
    executable_settings: {{
      events: {n_events},
      output: [ 
        {{ file: "{output_file}.hepevt", format: "hepevt", mode: "overwrite" }}
      ],
    }},
}}"""
    return config_text

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
    if len(sys.argv) != 3:
        print("Usage: python generate_snb_batch.py <job_id> <n_bursts>")
        sys.exit(1)
    
    job_id = int(sys.argv[1])
    n_bursts = int(sys.argv[2])
    n_events_per_burst = 200
    
    # Setup paths - use absolute path
    data_dir = Path("/afs/cern.ch/work/e/evilla/private/dune/marley-gen/data/snb_1000")
    data_dir.mkdir(parents=True, exist_ok=True)
    
    config_dir = Path("configs_temp")
    config_dir.mkdir(parents=True, exist_ok=True)
    
    # Storage for results
    all_directions = []
    all_electrons = []
    
    print(f"Job {job_id}: Generating {n_bursts} SNB events with {n_events_per_burst} ES interactions each")
    
    for burst_idx in range(n_bursts):
        global_idx = job_id * n_bursts + burst_idx
        
        # Generate random direction for this burst
        direction = generate_random_direction()
        all_directions.append(direction)
        
        print(f"  Burst {global_idx}: direction = ({direction[0]:.3f}, {direction[1]:.3f}, {direction[2]:.3f})")
        
        # Create config
        output_file = str(data_dir / f"snb_{global_idx:04d}")
        config_file = config_dir / f"config_{job_id}_{burst_idx}.js"
        
        config_text = create_marley_config(direction, output_file, n_events_per_burst)
        hepevt_file = output_file + ".hepevt"
        
        with open(config_file, 'w') as f:
            f.write(config_text)
        
        # Run MARLEY
        try:
            marley_exe = os.environ.get('MARLEY', '/afs/cern.ch/work/e/evilla/private/dune/marley-1.2.0') + '/build/marley'
            result = subprocess.run(
                [marley_exe, str(config_file)],
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                print(f"    ERROR: MARLEY failed for burst {global_idx}")
                print(result.stderr)
                all_electrons.append(None)
                continue
            
            # Extract electrons
            electrons = extract_electrons_from_hepevt(hepevt_file)
            
            if electrons is not None and len(electrons) > 0:
                all_electrons.append(electrons)
                print(f"    ✓ Extracted {len(electrons)} electrons")
                # Clean up HEPEVT file to save space
                os.remove(hepevt_file)
            else:
                print(f"    WARNING: No electrons extracted for burst {global_idx}")
                all_electrons.append(None)
            
        except subprocess.TimeoutExpired:
            print(f"    ERROR: MARLEY timeout for burst {global_idx}")
            all_electrons.append(None)
        except Exception as e:
            print(f"    ERROR: {e}")
            all_electrons.append(None)
    
    # Save results for this job
    output_file = data_dir / f"job_{job_id:03d}_results.npz"
    
    np.savez(
        output_file,
        directions=np.array(all_directions),
        electrons=np.array(all_electrons, dtype=object),
        job_id=job_id,
        n_bursts=n_bursts
    )
    
    print(f"\nJob {job_id}: Complete! Saved to {output_file}")
    print(f"  Total bursts: {n_bursts}")
    print(f"  Successful: {sum(1 for e in all_electrons if e is not None)}")
    
    # Clean up temp configs
    for f in config_dir.glob(f"config_{job_id}_*.js"):
        f.unlink()

if __name__ == "__main__":
    main()
