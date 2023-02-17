from md_analysis import trj
from md_analysis import atoms

cell= [
    [21.2944000000,        0.0000000000,        0.0000000000],
    [-4.6030371123,       20.7909480472,        0.0000000000],
    [-0.9719093466,       -1.2106211379,       15.1054299403]
]
# traj = trj.Trajectory(timestep = 5, data = "cage1-500K.xyz", index = ":", cell = cell)
traj = trj.Trajectory(timestep = 5, data = "cage1-500K.xyz", index = slice(1000,11001,1), cell = cell)

from recordclass import make_dataclass, dataobject, astuple, asdict
water = make_dataclass("Water", 'O H1 H2')
waters=[[14,15,16,], # Away from cage
        [17,18,19,], # Hydrogen bonded to amine
        [143,144,145,], # Hydrogen bonded to phenol
        [146,147,148,], # Showed significant diffusion
        [272,273,274,], # Hydrogen bonded to phenol
        [275,276,277,], # Showed significant diffusion
        [401,402,403,], # Away from cage
        [404,405,406,], # Showed significant diffusion
        ]
for i in range(len(waters)):
    globals()[f"water{i+1}"] = water(*waters[i])
    
TFA = make_dataclass("TFA", "F1 F2 F3 O1 O2 C1 C2")
TFAs = [
    [0,1,2,3,4,5,6,], # Free
    [7,8,9,10,11,12,13,], # Bound
    [129,130,131,132,133,134,135,], # Bound
    [136,137,138,139,140,141,142,], # Bound
    [258,259,260,261,262,263,264,], # Bound
    [265,266,267,268,269,270,271,], # Bound
    [387,388,389,390,391,392,393,], # Free
    [394,395,396,397,398,399,400,], # Bound
]
for i in range(len(TFAs)):
    globals()[f"TFA{i+1}"] = TFA(*TFAs[i])
    
TFA_F_O = make_dataclass("TFA_F_O", 'F1 F2 F3 O1 O2')
temp_TFA_F_O = TFA_F_O(*list(range(5)))
for i in range(len(TFAs)):
    keys = list(asdict(temp_TFA_F_O).keys())
    globals()[f"TFA_F_O{i+1}"] = TFA_F_O(*[asdict(globals()[f"TFA{i+1}"])[key] for key in keys])
    
phenol = make_dataclass("Phenol", "O H")
phenols = [
    [24, 25], #
    [411, 412], #
    [153, 154], # H-bonded to water O 272
    [282, 283], # H-bonded to water O 143
]
      
for i in range(len(phenols)):
    globals()[f"phenol{i+1}"] = phenol(*phenols[i])
    
amine_2H = make_dataclass("Amine", "N H1 H2")
amine_1H = make_dataclass("Amine", "N H")
amines = [
    [410, 284, 285], #H1 H-bond to TFA O 397
    [409, 128], # H-bond to TFA O 261
    [407, 28, 29], #H1 H-bond to TFA O 261, H2 H-bonded to TFA O 269
    [152, 26, 27], #
    [149, 286, 287], #
    [151, 386,], #
    [280, 257], #
    [278, 157, 158], #
    [281, 413, 414], #
    [23, 155, 156], # H1 155 H-bonded to TFA O 10 H2 156 H-bonded to water O 17
    [22, 515], # H 149 H-bonded to TFA O 132
    [20, 415, 416], # H1 415 H-bonded to TFA O 132 H2 416 H-bonded to TFA O 140
]

for i, amine  in enumerate(amines):
    if len(amine) ==2:
        globals()[f"amine{i+1}"] = amine_1H(*amines[i])
    else:
        globals()[f"amine{i+1}"] = amine_2H(*amines[i])
        
atoms_dict = [phenol(*phenols[i]) for i in range(len(phenols))] + \
    [amine_1H(*amines[i]) if len(amine) ==2 else amine_2H(*amines[i]) for i, amine  in enumerate(amines)] + \
    [TFA_F_O(*[asdict(globals()[f"TFA{i+1}"])[key] for key in list(asdict(temp_TFA_F_O).keys())]) for i in range(len(TFAs))] + \
    [water(*waters[i]) for i in range(len(waters))]

def retrieve_symbol(string):
    """function to remove numbers in a string, so that the atom dict keys can be converted to chemical symbols"""
    return ''.join([i for i in string if not i.isdigit()])

def find_d_a(atoms_dict):
    """
    Creat donor, acceptor dicts from atoms_dicts. Include one water to creat donor, acceptor for specified water.
    """
    donors = []
    acceptors = []
    nd = make_dataclass("N_donor", "N H")
    od = make_dataclass("O_donor", "O H")
    oa = make_dataclass("O_acceptor", "O")
    for a_dict in atoms_dict:
        symbols = [retrieve_symbol(key) for key in list(asdict(a_dict).keys())] # retrieve chemical symbol from atoms_dict keys
        if symbols[0] == "F":
            acceptors += (a_dict,)
        elif symbols[0] == "N":
            if len(symbols) == 2:
                donors += (a_dict,)
            else:
                combs = [[asdict(a_dict)["N"], asdict(a_dict)[x]] for x in list(asdict(a_dict).keys()) if x != "N"]
                donors += (nd(*combs[0]),)
                donors += (nd(*combs[1]),)

        elif symbols[0] == "O":
            if symbols[0] == symbols[1]:
                acceptors += (a_dict,)
            elif len(symbols) == 2 and symbols[0] != symbols[1]:
                donors += (a_dict,)
            else:
                combs = [[asdict(a_dict)["O"], asdict(a_dict)[x]] for x in list(asdict(a_dict).keys()) if x != "O"]
                donors += (od(*combs[0]),)
                donors += (od(*combs[1]),)
                acceptors += (oa(asdict(a_dict)["O"]),)
    print(f"Done! {len(donors)} donors and {len(donors)} acceptors were found!\n" + "Acceptors:\n", acceptors, "\n Donors:\n", donors)
    return donors, acceptors

donors, acceptors = find_d_a(atoms_dict = atoms_dict)

def pair_d_a(donors, acceptors):
    """
    Make d_a_pairs and d_h_a_pairs from the donors and acceptors dict.
    """
    d_a_pairs = []
    d_h_a_pairs = []
    for donor in donors:
        d_atom = [{x: asdict(donor)[x]} for x in list(asdict(donor).keys()) if not x.startswith("H")]
        for acceptor in acceptors:
            keys = list(asdict(acceptor).keys())
            if len(keys) == 1:
                symbol = list(d_atom[0].keys())[0]
                if symbol == "O" and d_atom[0][symbol] == asdict(acceptor)[keys[0]]: # exclude combinations where donor and acceptor are the same oxygen atom.
                    print(f"Donor and acceptor is the same atom: {symbol} {asdict(acceptor)[keys[0]]}, skipped.")
                    continue
                else:
                    d_a_pairs += [[d_atom[0], asdict(acceptor)],]
                    d_h_a_pairs += [[asdict(donor), asdict(acceptor)],]
            else:
                for key in keys:
                    d_a_pairs += [[d_atom[0], {key: asdict(acceptor)[key]}],]
                    d_h_a_pairs += [[asdict(donor), {key: asdict(acceptor)[key]}],]
    print(f"Unique donoar-acceptor/donoar-hydrogen-acceptor combinations: {len(d_a_pairs)}")
    print(f"Sample donor-acceptor pair: {d_a_pairs[0]}\nSample donor-hydrogen-acceptor pair: {d_h_a_pairs[0]}")
    return d_a_pairs, d_h_a_pairs

d_a_pairs, d_h_a_pairs = pair_d_a(donors = donors, acceptors = acceptors)

def at_dict2pos(at_dict, frame):
    """
    Converts the list of at_dict to list of positions
    """
    for item in at_dict:
        indices = (list(item.values()))
    return [frame.pos[i] for i in indices]

def distance(atom_pair, frame):
    """
    Calculates the distance of a atom pair.
    """
    indices = [list(atom_dict.values()) for atom_dict in atom_pair]
    indices = [item for sublist in indices for item in sublist] # Flatten the list of lists to a list of integers
    # positions = [frame.pos[i] for i in indices]
    # # print(positions)
    # sqdst = np.square(positions[1] - positions[0]).sum()
    # dist = np.sqrt(sqdst)
    return frame.dist(indices[0], indices[1], mic = True)

from md_analysis.data import vdW_R
import numpy as np
def theta(vec1, vec2):
    """
    Calculate the angle between two vectors.
    """
    return np.arccos(vec1.dot(vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))*180/np.pi
def res_h(d_a_pairs, d_h_a_pairs, frames):
    """
    Calculates the H-bond information in each frame of the trajectory.
    """
    results = []
    for i, frame in enumerate(frames):
        hbonds = [] #list to store all hbonds in current frame
        for d_a_pair, d_h_a_pair in zip(d_a_pairs, d_h_a_pairs):
            # Sum of van der Waals radii
            symbols = [list(atom_dict.keys()) for atom_dict in d_a_pair]
            symbols = [item for sublist in symbols for item in sublist] # flatten the list of lists
            symbols = [retrieve_symbol(symbol) for symbol in symbols] # Remove numbers and get chemical symbol from atoms_dict to be ready to pass to the vdW radii dict
            vdW_sum = vdW_R[symbols[0]] + vdW_R[symbols[1]]
            
            # not a H-bond of D-A distance is greater than their vdW radii times 1.05, 1.05 to take bond length change during MD simulation.
            r_d_a = distance(d_a_pair, frame)  # calculate the D-A distance
            if r_d_a < 1.02 * vdW_sum:
                # calculate the D-H⋅⋅⋅A angle
                d_h_pos = at_dict2pos ([d_h_a_pair[0],], frame) # the positions of donor and hydrogen
                a_pos = at_dict2pos([d_h_a_pair[1],], frame) # the positions of acceptor
                d_h_vec = np.asarray(d_h_pos[1], dtype = float) - np.asarray(d_h_pos[0], dtype = float)
                a_h_vec = np.asarray(d_h_pos[1], dtype = float) - np.asarray(a_pos[0], dtype = float)
                d_h_a_ang = theta(d_h_vec, a_h_vec) # angle
                d_h = distance([d_h_a_pair[0],], frame) # calculate the D-H length
            
                # the D-H⋅⋅⋅A angle criteria used: the D-H⋅⋅⋅A angle is close to a right angle refer to the D-H⋅⋅⋅A angle - R(D⋅⋅⋅A) plot
                # an angle range is included considering the oscillation of bond lenghth and anlgle
                if d_h_a_ang >= (np.rad2deg(np.arctan2(r_d_a, d_h)) + 180)*3/8:
                # if d_h_a_ang >= 90:
                    # Store current H-bond
                    hbonds.append(
                          {
                              "donor": d_h_a_pair[0],
                              "acceptor": d_h_a_pair[1],
                              "R(D-A)": r_d_a,
                              "DHA_ang": d_h_a_ang,
                              "r(D-H)": d_h,
                          }
                    )
                
        results.append({f"frame": i, "n_hbonds": len(hbonds), "hbonds": hbonds})
    return results

print("Resolving H-bonds... ...")
hbonds = res_h(d_a_pairs=d_a_pairs, d_h_a_pairs=d_h_a_pairs, frames = traj)
import json
with open("hbonds-cage1-500K.json", "w") as json_file:
    json.dump(hbonds, json_file)
json_file.close()
