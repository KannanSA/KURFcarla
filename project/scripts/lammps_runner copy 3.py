# scripts/lammps_runner.py
import os
import pandas as pd
from ase import Atoms
from ase.io import read, write

# Path to GAP-20 XML file
GAP_XML = '/Users/kannansekarannuradha/Documents/Projects/KURFcarla/project/results/Carbon_GAP_20.xml'  # Use absolute path

# Reference energy per atom (eV) for graphite
E_REF = -7.37

def run_binding_energy(xyz_file, out_dir, lammps_executable='lammps'):
    """
    Runs a single-point energy calculation using the GAP-20 potential
    on the cluster defined in xyz_file, computes binding energy per atom,
    writes log, and returns results dict.
    """
    atoms = read(xyz_file)
    n_atoms = len(atoms)
    base = os.path.splitext(os.path.basename(xyz_file))[0]
    os.makedirs(out_dir, exist_ok=True)

    # Write LAMMPS data file
    data_file = os.path.join(out_dir, base + '.data')
    write(data_file, atoms, format='lammps-data', units='real', atom_style='atomic')

    # Select GAP-20 force field
    params = {
        'pair_style': 'quip',
        'pair_coeff': [f'* * {GAP_XML} "Potential" C'],
        'mass': ['1 12.011']
    }
    from ase.calculators.lammpsrun import LAMMPS
    try:
        calc = LAMMPS(**params, files=[GAP_XML], command=lammps_executable)
        atoms.calc = calc
        E_cluster = atoms.get_potential_energy()
        BE_per_atom = (n_atoms * E_REF - E_cluster) / n_atoms
    except RuntimeError as e:
        msg = str(e)
        if (
            "pair style 'quip'" in msg
            or "ML-QUIP package" in msg
            or "Unrecognized pair style 'quip'" in msg
            or "LAMMPS exited" in msg
        ):
            print("ERROR: Your LAMMPS binary does not have the ML-QUIP package enabled, which is required for 'pair_style quip'.")
            print("Please install or compile a LAMMPS binary with the ML-QUIP package enabled (see https://github.com/libAtoms/lammps-quip).")
            print("Original error message:\n", msg)
            raise
        elif "Failed to retrieve any thermo_style-output" in msg:
            print("ERROR: LAMMPS did not produce any thermo_style output. This usually means LAMMPS failed to run properly, likely due to a missing or misconfigured package (such as ML-QUIP).")
            print("Check your LAMMPS build and input files.")
            print("Original error message:\n", msg)
            raise
        else:
            raise

    # Write log
    log_file = os.path.join(out_dir, base + '.log')
    with open(log_file, 'w') as f:
        f.write(f"TotalEnergy {E_cluster}\n")
        f.write(f"BE_per_atom {BE_per_atom}\n")

    return {
        'cluster': base,
        'n_atoms': n_atoms,
        'E_cluster': E_cluster,
        'BE_per_atom': BE_per_atom
    }

def batch_run(cluster_dir, out_dir, results_csv, lammps_executable='lammps'):
    """
    Processes all .xyz files in cluster_dir, runs binding-energy calc,
    and saves results to results_csv.
    """
    results = []
    for fname in os.listdir(cluster_dir):
        if not fname.endswith('.xyz'):
            continue
        xyz_path = os.path.join(cluster_dir, fname)
        print(f"Running on {fname}...")
        res = run_binding_energy(xyz_path, out_dir, lammps_executable)
        results.append(res)

    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(results_csv), exist_ok=True)
    df.to_csv(results_csv, index=False)
    print(f"Results saved to {results_csv}")
    return df

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run GAP-20 binding energy calculations via ASE/LAMMPS')
    # parser.add_argument('--xyz', type=str, help='Path to a single XYZ file')
    parser.add_argument('--xyz', type=str, default='/Users/kannansekarannuradha/Documents/Projects/KURFcarla/project/data/clusters/C-2158988-7926-291.xyz', help='project/data/clusters/C-2158988-7926-291.xyz')
    parser.add_argument('--dir', type=str, default='project/data/clusters', help='Directory of XYZ cluster files')
    parser.add_argument('--out', type=str, default='project/out', help='Output directory for data and logs')
    parser.add_argument('--results', type=str, default='/Users/kannansekarannuradha/Documents/Projects/KURFcarla/project/results/Carbon_GAP_20.xml', help='CSV path to save aggregated results')
    parser.add_argument('--lammps', type=str, default='/Users/kannansekarannuradha/Documents/Projects/KURFcarla/lammps/src/lmp_mpi', help='Path to the LAMMPS executable or directory containing it')
    # parser.add_argument('--lammps', type=str, default='/Users/kannansekarannuradha/opt/anaconda3/envs/lammps-ase/bin/lmp', help='Path to the LAMMPS executable or directory containing it')
#/usr/local/Cellar/lammps/20240829-update2/bin/lmp_mpi
#/usr/local/Cellar/lammps/20240829-update2/bin/lmp_serial
    args = parser.parse_args()

    lmp_executable = args.lammps
    if os.path.isdir(lmp_executable):
        lmp_executable = os.path.join(lmp_executable, 'lmp')
    # os.environ["LAMMPS_COMMAND"] = lmp_executable  # <-- REMOVE THIS LINE

    # from ase.calculators.lammpsrun import LAMMPS  # <-- REMOVE THIS LINE

    # Check if the LAMMPS executable exists
    if not os.path.exists(lmp_executable):
        raise FileNotFoundError(f"LAMMPS executable not found at {lmp_executable}. Check your --lammps argument.")
    
    if args.xyz:
        os.makedirs(os.path.dirname(args.results), exist_ok=True)
        try:
            res = run_binding_energy(args.xyz, args.out, lmp_executable)
            df = pd.DataFrame([res])
            df.to_csv(args.results, index=False)
            print(f"Result for {args.xyz}:\n", df)
        except RuntimeError as e:
            msg = str(e)
            if (
                "pair style 'quip'" in msg
                or "ML-QUIP package" in msg
                or "Unrecognized pair style 'quip'" in msg
                or "LAMMPS exited" in msg
            ):
                print("ERROR: Your LAMMPS binary does not have the ML-QUIP package enabled, which is required for 'pair_style quip'.")
                print("Please install or compile a LAMMPS binary with the ML-QUIP package enabled (see https://github.com/libAtoms/lammps-quip).")
                print("Original error message:\n", msg)
            else:
                raise
    else:
        try:
            df = batch_run(args.dir, args.out, args.results, lmp_executable)
            print(df)
        except RuntimeError as e:
            msg = str(e)
            if (
                "pair style 'quip'" in msg
                or "ML-QUIP package" in msg
                or "Unrecognized pair style 'quip'" in msg
                or "LAMMPS exited" in msg
            ):
                print("ERROR: Your LAMMPS binary does not have the ML-QUIP package enabled, which is required for 'pair_style quip'.")
                print("Please install or compile a LAMMPS binary with the ML-QUIP package enabled (see https://github.com/libAtoms/lammps-quip).")
                print("Original error message:\n", msg)
            else:
                raise