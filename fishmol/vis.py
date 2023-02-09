from ase import Atoms as ase_Atoms
from ase.visualize import view as ase_view
from IPython.display import HTML, display
from tempfile import NamedTemporaryFile

def view(atoms, viewer = "ngl"):
    atoms = ase_Atoms(symbols = atoms.symbs, positions = atoms.pos, cell = atoms.cell.lattice)
    ase_view(atoms)

def render_atoms(atoms):
    """
    render the ase atoms object.
    A: centre view
    E: rotate
    Double click on atom to bring the target atom into centre of view
    Roll to zoom in/out
    """
    atoms = ase_Atoms(symbols = atoms.symbs, positions =  atoms.pos)
    with NamedTemporaryFile('r+', suffix='.html') as ntf:
        atoms.write(ntf.name, format='html')
        ntf.seek(0)
        html = ntf.read()
        display(HTML(html))
    return None