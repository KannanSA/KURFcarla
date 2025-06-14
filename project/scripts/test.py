from ase import Atoms
from ase.calculators.quip import Quip

atoms = Atoms('Si2', positions=[[0,0,0], [0,0,2.35]])
calc = Quip(potential='IP SW')
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
print(f'Energy: {energy} eV')
