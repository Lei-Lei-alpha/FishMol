import itertools
from recordclass import make_dataclass
from fishmol.data import val_R

def cluster(atoms, mic = False):
    """
    Clustering the atoms object by calculating the bonded atoms. Return a list of Molecule dataclass that the formula and list of at_idx can be extracted by .formula and .at_idx methods.
    """
    # Make a dataclass to store the info of clusted molecules
    molecule = make_dataclass("Molecule", "formula at_idx")
    # All unique combinations of indices for atoms 
    combinations = list((i, j) for ((i,_),(j,_)) in itertools.combinations(enumerate(list(atoms.symbs)), 2))
    bonded = []
    for comb in combinations:
        if mic:
            dist = atoms.dist(comb[0], comb[1], mic = True)
        else:
            dist = atoms.dist(comb[0], comb[1], mic = False)
        # The atom pairs are bonded if their distance is less than 1.05 times their valency radii
        if dist < 1.05 * (val_R[atoms.symbs[comb[0]]] + val_R[atoms.symbs[comb[1]]]):
            bonded.append(comb)
    
    # Cluster the bonded atom pairs
    i = 0
    while True:
        globals()[f"mol{i}"] = []

        if globals()[f"mol{i}"] == []:
            if len(bonded) == 1:
                globals()[f"mol{i}"] = bonded
                bonded = []
            else:
                a = [bonded[x] for x in range(len(bonded)) if any(ele in bonded[x] for ele in bonded[0])]
                bonded = [ele for ele in bonded if not ele in a]
                globals()[f"mol{i}"] += a

        while any(ele in itertools.chain.from_iterable(bonded) for ele in itertools.chain.from_iterable(globals()[f"mol{i}"])):
            # print(share_node)
            a = [bonded[x] for x in range(len(bonded)) if any(ele in bonded[x] for ele in itertools.chain.from_iterable(globals()[f"mol{i}"]))]
            bonded = [ele for ele in bonded if not ele in a]
            globals()[f"mol{i}"] += a

        if bonded !=[]:
            i += 1
            continue

        if bonded == []:
            mols = [list(set(itertools.chain.from_iterable(globals()[f"mol{x}"]))) for x in range(i+1)]
            break
    # Resolve the formula by counting the number of atoms in each molecule
    symbols = [atoms.symbs[mol].tolist() for mol in mols]
    s_list = [list(set(symb)) for symb in symbols] # list of symbols without duplicates
    counts = [[str(a.count(n)) for n in b] for a,b in zip(symbols, s_list)] # count number of atoms for each symbol
    symb_num_comb = [list(itertools.chain.from_iterable(zip(a, b))) for a, b in zip(s_list, counts)]
    # symb_num_comb = [list(filter(('1').__ne__, x)) for x in symb_num_comb] # Drop '1' from chemical formula
    formula = ["".join(x) for x in symb_num_comb]
    mols = [molecule(a, b) for a, b in zip(formula, mols)]
    return mols
