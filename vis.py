from ase import Atoms as ase_Atoms
from ase.visualize import view as ase_view

def view(atoms, viewer = "ngl"):
    atoms = ase_Atoms(symbols = atoms.symbs, positions = atoms.pos, cell = atoms.cell.lattice)
    ase_view(atoms)