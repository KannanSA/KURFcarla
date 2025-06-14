# scripts/projector.py
import os
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorboard.plugins import projector
from ase.io import read


def generate_embeddings(results_csv, clusters_dir, logdir):
    """
    Reads results_csv, computes radius of gyration for each cluster,
    builds an embedding of [n_atoms, BE_per_atom, Rg], writes metadata,
    and configures TensorBoard projector in logdir.
    """
    os.makedirs(logdir, exist_ok=True)

    # Load results
df = pd.read_csv(results_csv)

    # Compute radius of gyration
def compute_rg(fname):
    atoms = read(os.path.join(clusters_dir, fname + '.xyz'))
    com = atoms.get_center_of_mass()
    coords = atoms.get_positions()
    return np.sqrt(np.mean(np.sum((coords - com)**2, axis=1)))

rgs = [compute_rg(c) for c in df['cluster']]
df['Rg'] = rgs

# Prepare embeddings array
embeddings = np.vstack([df['n_atoms'], df['BE_per_atom'], df['Rg']]).T

# Save metadata.tsv
metadata_path = os.path.join(logdir, 'metadata.tsv')
with open(metadata_path, 'w') as f:
    f.write('cluster	n_atoms	BE_per_atom	Rg
')
    for _, row in df.iterrows():
        f.write(f"{row['cluster']}	{row['n_atoms']}	{row['BE_per_atom']:.6f}	{row['Rg']:.6f}
")

# Create TensorFlow variable
embedding_var = tf.Variable(embeddings, name='cluster_embeddings')

# Save checkpoint
checkpoint = tf.train.Checkpoint(embeddings=embedding_var)
checkpoint_prefix = os.path.join(logdir, 'embeddings.ckpt')
checkpoint.save(checkpoint_prefix)

# Configure projector
config = projector.ProjectorConfig()
emb = config.embeddings.add()
emb.tensor_name = embedding_var.name
emb.metadata_path = 'metadata.tsv'
projector.visualize_embeddings(os.path.abspath(logdir), config)

if __name__ == '__main__':
    generate_embeddings(
        'project/results/binding_energies.csv',
        'project/data/clusters',
        'project/results/projector'
    )
