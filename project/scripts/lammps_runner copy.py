import os
import subprocess
from ase import Atoms
from ase.io import read, write
from ase.calculators.lammpsrun import LAMMPS  # changed import to use lammpsrun

# 1) Read cluster structure, write LAMMPS data
def load_cluster(path):
    return read(path)

def write_lammps_data(atoms, out_data):
    write(out_data, atoms, format='lammps-data', units='real',
          atom_style='atomic')

# 2) Run LAMMPS via ASE LAMMPS calculator
def run_lammps(data_file, parameters, logfile):
    """Runs LAMMPS via ASE calculator. 
    parameters: dict of LAMMPS parameters (pair_style, coeffs, etc.)"""
    calc = LAMMPS(parameters=parameters, files=['in.template'])
    atoms = read(data_file, format='lammps-data')
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()  # eV
    with open(logfile, 'w') as f:
        f.write(f"TotalEnergy {energy}\n")
    return energy

# 3) Batch over directory
def batch_run(cluster_dir, out_dir, lammps_params, lammps_exec='lmp'):
    os.makedirs(out_dir, exist_ok=True)
    results = []
    for fname in os.listdir(cluster_dir):
        if not fname.endswith(('.xyz','.pdb')): continue
        base = os.path.splitext(fname)[0]
        struct = os.path.join(cluster_dir, fname)
        data = os.path.join(out_dir, base+'.data')
        log  = os.path.join(out_dir, base+'.log')
        atoms = load_cluster(struct)
        write_lammps_data(atoms, data)
        E = run_lammps(data, lammps_params, log)
        results.append({'cluster': base, 'n_atoms': len(atoms), 'E_cluster': E})
    return results

if __name__ == '__main__':
    # Example parameters: Tersoff for carbon
    params = {
      'pair_style': 'tersoff',
      'pair_coeff': ['* * SiC.tersoff C']
    }
    out = batch_run('project/data/clusters', 'project/out', params)
    # Save raw energies
    import pandas as pd
    df = pd.DataFrame(out)
    df.to_csv('project/results/results.csv', index=False)  # changed output to a file
