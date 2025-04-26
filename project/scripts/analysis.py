import pandas as pd
import matplotlib.pyplot as plt

# 1) Load raw cluster energies
df = pd.read_csv('results/raw_energies.csv')

# 2) Reference: energy per C atom in bulk (from separate sim)
E_ref = -7.37  # [eV/atom] for graphite

# 3) Compute binding energy per atom:  
# BE = (n * E_ref – E_cluster)/n
df['BE_per_atom'] = (df['n_atoms']*E_ref - df['E_cluster'])/df['n_atoms']

# 4) Save results
df.to_csv('results/binding_energies.csv', index=False)

# 5) Plot trend vs size
plt.figure()
plt.scatter(df['n_atoms'], df['BE_per_atom'])
plt.xlabel('Number of atoms')
plt.ylabel('Binding energy per atom (eV)')
plt.title('Binding energy trend for carbon nanoclusters')
plt.tight_layout()
plt.savefig('results/binding_energy_vs_size.png', dpi=300)
plt.show()
