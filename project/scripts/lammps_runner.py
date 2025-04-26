# scripts/lammps_runner.py
import os
import pandas as pd
from ase import Atoms
from ase.io import read, write
from ase.calculators.lammpsrun import LAMMPS

# Path to your GAP-20 XML potential file
GAP_XML = 'project/results/Carbon_GAP_20 copy.xml'
# GAP_XML = 'project/results/Carbon_GAP_20.xml'


def batch_run(cluster_dir, out_dir):
    """
    Read all .xyz/.pdb files in cluster_dir, convert to LAMMPS data,
    run single-point energy calc using GAP-20 via LAMMPS, save logs.
    Returns DataFrame with cluster name, atom count, and energy.
    """
    os.makedirs(out_dir, exist_ok=True)
    results = []
    files = [f for f in os.listdir(cluster_dir) if f.endswith(('.xyz', '.pdb'))]
    if not files:
        print(f"Warning: No cluster files found in {cluster_dir}.")
        return pd.DataFrame(results)

    for fname in files:
        base = os.path.splitext(fname)[0]
        data_file = os.path.join(out_dir, base + '.data')
        log_file  = os.path.join(out_dir, base + '.log')

        # Read cluster and write LAMMPS-format data
        atoms = read(os.path.join(cluster_dir, fname))
        write(data_file, atoms, format='lammps-data', units='real', atom_style='atomic')

        # Setup QUIP/GAP20 with LAMMPS
        params = {
            'pair_style': 'quip',
            'pair_coeff': [f'* * {GAP_XML} C'],
            'mass':        ['1 12.011'],
        }
        calc = LAMMPS(parameters=params, files=[GAP_XML])
        atoms.set_calculator(calc)

        # Compute energy
        try:
            E = atoms.get_potential_energy()
        except Exception as err:
            print(f"Error computing energy for {fname}: {err}")
            E = float('nan')

        # Write log
        with open(log_file, 'w') as f:
            f.write(f"TotalEnergy {E}\n")  # fixed: entire string on one line

        results.append({'cluster': base,
                        'n_atoms': len(atoms),
                        'E_cluster': E})

    df = pd.DataFrame(results)
    print(df)
    return df

if __name__ == '__main__':
    clusters_dir = 'project/data/clusters'
    out_dir      = 'project/out'
    results_df   = batch_run(clusters_dir, out_dir)
    # Save final results
    os.makedirs('project/results', exist_ok=True)
    results_df.to_csv('project/results/results.csv', index=False)

# scripts/analysis.py
import pandas as pd
import matplotlib.pyplot as plt

# Load energies from the correct CSV file
df = pd.read_csv('project/results/results.csv')  # updated file path
# Reference energy per atom (eV) for graphite
E_ref = -7.37
# Compute binding energy per atom
df['BE_per_atom'] = (df['n_atoms'] * E_ref - df['E_cluster']) / df['n_atoms']

# Save binding energies
df.to_csv('results/binding_energies.csv', index=False)

# Plot trend vs cluster size
plt.figure()
plt.scatter(df['n_atoms'], df['BE_per_atom'])
plt.xlabel('Number of atoms')
plt.ylabel('Binding energy per atom (eV)')
plt.title('Binding energy trend (GAP-20 carbon)')
plt.tight_layout()
plt.savefig('results/binding_energy_vs_size.png', dpi=300)
plt.show()


# scripts/ovito_viz.py
import os
from ovito.io import import_file, export_file  # added import for ovito

# Export a snapshot image for each LAMMPS dump
def export_snapshot(dump_file, out_image, frame=0, dpi=300, size=(800,600)):
    pipeline = import_file(dump_file)
    export_file(pipeline, out_image, 'png', frame=frame, size=size, dpi=dpi)

if __name__ == '__main__':
    os.makedirs('results/ovito', exist_ok=True)
    for fname in os.listdir('out/'):
        if fname.endswith('.dump'):
            in_dump = os.path.join('out/', fname)
            out_img  = os.path.join('results/ovito', fname.replace('.dump', '.png'))
            export_snapshot(in_dump, out_img)