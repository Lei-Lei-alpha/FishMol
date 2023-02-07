import mmap
import numpy as np

# van der Waals radii
vdW_R = {"O": 1.52, "N": 1.55, "F": 1.47,}

def read_traj(traj, prop = None, ):
    """
    read the xyz trajctory file, and remove the headers: line1: number of atoms in the system, line2: Properties of the system, pbc, simu step, time, a,b,c etc.
    prop: input the start of second header if it is not "Properties"
    """
    with open(traj) as f:
        mm = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
        contents = [line.rstrip().decode('utf8') for line in iter(mm.readline, b"")]
        mm.close()
    f.close()
    header = contents[0]
    natoms = int(header)
    nframes = contents.count(contents[0])
    if prop == None:
        prop = contents[1].split("=")
    contents = [line for line in contents if not (line.startswith(header) or line.startswith(prop[0]))]
    contents = [contents[x: x + natoms] for x in range(0, len(contents), natoms)]
    return natoms, nframes, contents

def find_d_a(atoms_dict):
    """
    Creat donor, acceptor dicts from atoms_dicts. Include one water to creat donor, acceptor for specified water.
    """
    donors = []
    acceptors = []
    for a_dict in atoms_dict:
        symbols = [retrieve_symbol(key) for key in list(a_dict.keys())] # retrieve chemical symbol from atoms_dict keys
        if symbols[0] == "F":
            acceptors += (a_dict,)
        elif symbols[0] == "N":
            if len(symbols) == 2:
                donors += (a_dict,)
            else:
                combs = [{"N": a_dict["N"], x: a_dict[x]} for x in list(a_dict.keys()) if x != "N"]
                donors += (combs[0], combs[1])

        elif symbols[0] == "O":
            if symbols[0] == symbols[1]:
                acceptors += (a_dict,)
            elif len(symbols) == 2 and symbols[0] != symbols[1]:
                donors += (a_dict,)
            else:
                combs = [{"O": a_dict["O"], x: a_dict[x]} for x in list(a_dict.keys()) if x != "O"]
                donors += (combs[0], combs[1])
                acceptors += ({"O": a_dict["O"]},)
    print(f"Done! {len(donors)} donors and {len(donors)} acceptors were found!\n" + "Acceptors:\n", acceptors, "\n Donors:\n", donors)
    return donors, acceptors

def pair_d_a(donors, acceptors):
    """
    Make d_a_pairs and d_h_a_pairs from the donors and acceptors dict.
    """
    d_a_pairs = []
    d_h_a_pairs = []
    for donor in donors:
        d_atom = [{x: donor[x]} for x in list(donor.keys()) if not x.startswith("H")]
        for acceptor in acceptors:
            keys = list(acceptor.keys())

            if len(keys) == 1:
                symbol = list(d_atom[0].keys())[0]
                if symbol == "O" and d_atom[0][symbol] == acceptor[keys[0]]: # exclude combinations where donor and acceptor are the same oxygen atom.
                    print(f"Donor and acceptor is the same atom: {symbol} {acceptor[keys[0]]}, skipped.")
                    continue
                else:
                    d_a_pairs += [[d_atom[0], acceptor],]
                    d_h_a_pairs += [[donor, acceptor],]
            else:
                for key in keys:
                    d_a_pairs += [[d_atom[0], {key: acceptor[key]}],]
                    d_h_a_pairs += [[donor, {key: acceptor[key]}],]
    print(f"Unique donoar-acceptor/donoar-hydrogen-acceptor combinations: {len(d_a_pairs)}")
    print(f"Sample donor-acceptor pair: {d_a_pairs[0]}\nSample donor-hydrogen-acceptor pair: {d_h_a_pairs[0]}")
    return d_a_pairs, d_h_a_pairs

def theta(vec1, vec2):
    """
    Calculate the angle between two vectors.
    """
    return np.arccos(vec1.dot(vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))*180/np.pi

def retrieve_symbol(string):
    """function to remove numbers in a string, so that the atom dict keys can be converted to chemical symbols"""
    return ''.join([i for i in string if not i.isdigit()])

def at_dict2pos(at_dict, frame):
    """
    Converts the list of at_dict to list of positions
    """
    for item in at_dict:
        indices = (list(item.values()))
        positions = [frame[x].split(" ")[-3:] for x in indices]
    return positions

def distance(atom_pair, frame):
    """
    Calculates the distance of a atom pair.
    """
    indices = [list(atom_dict.values()) for atom_dict in atom_pair]
    indices = [item for sublist in indices for item in sublist] # Flatten the list of lists to a list of integers
    positions = [frame[x].split(" ")[-3:] for x in indices]
    if len(positions[0]) != 3 or len(positions[1]) != 3:
        print("The dimension of position is not correct in the current positions:\n", positions)

    sqdst = np.square(np.array(positions[1], dtype = float) - np.array(positions[0], dtype = float)).sum()
    dist = np.sqrt(sqdst)
    return dist

def res_h(d_a_pairs, d_h_a_pairs, frames, water_d = None, water_a = None):
    """
    Calculates the H-bond information in each frame of the trajectory.
    """
    results = []
    if water_d != None:
        for i, frame in enumerate(frames):
            hbonds = [] #list to store all hbonds in current frame
            for d_a_pair, d_h_a_pair in zip(d_a_pairs, d_h_a_pairs):
                # Keep water participated H-bonds only
                if d_h_a_pair[0] in water_d or d_h_a_pair[1] in water_a:
                    # Sum of van der Waals radii
                    symbols = [list(item.keys()) for item in d_a_pair]
                    symbols = [item for sublist in symbols for item in sublist] # flatten the list of lists
                    symbols = [retrieve_symbol(symbol) for symbol in symbols] # Remove numbers and get chemical symbol from atoms_dict to be ready to pass to the vdW radii dict
                    vdW_sum = vdW_R[symbols[0]] + vdW_R[symbols[1]]

                    # not a H-bond of D-A distance is greater than their vdW radii times 1.05, 1.05 to take bond length change during MD simulation.
                    r_d_a = distance(d_a_pair, frame)  # calculate the D-A distance
                    if r_d_a <= 1.02 * vdW_sum:
                        # calculate the D-H⋅⋅⋅A angle
                        d_h_pos = at_dict2pos ([d_h_a_pair[0],], frame) # the positions of donor and hydrogen
                        a_pos = at_dict2pos([d_h_a_pair[1],], frame) # the positions of acceptor
                        d_h_vec = np.array(d_h_pos[1], dtype = float) - np.array(d_h_pos[0], dtype = float)
                        a_h_vec = np.array(d_h_pos[1], dtype = float) - np.array(a_pos[0], dtype = float)
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
            
    else:
        for i, frame in enumerate(frames):
            hbonds = [] #list to store all hbonds in current frame
            for d_a_pair, d_h_a_pair in zip(d_a_pairs, d_h_a_pairs):
                # Sum of van der Waals radii
                symbols = [list(item.keys()) for item in d_a_pair]
                symbols = [item for sublist in symbols for item in sublist] # flatten the list of lists
                symbols = [retrieve_symbol(symbol) for symbol in symbols] # Remove numbers and get chemical symbol from atoms_dict to be ready to pass to the vdW radii dict
                vdW_sum = vdW_R[symbols[0]] + vdW_R[symbols[1]]

                # not a H-bond of D-A distance is greater than their vdW radii times 1.02, 1.02 to take bond length change during MD simulation.
                r_d_a = distance(d_a_pair, frame)  # calculate the D-A distance
                if r_d_a <= 1.02 * vdW_sum:
                    # calculate the D-H⋅⋅⋅A angle
                    d_h_pos = at_dict2pos ([d_h_a_pair[0],], frame) # the positions of donor and hydrogen
                    a_pos = at_dict2pos([d_h_a_pair[1],], frame) # the positions of acceptor
                    d_h_vec = np.array(d_h_pos[1], dtype = float) - np.array(d_h_pos[0], dtype = float)
                    a_h_vec = np.array(d_h_pos[1], dtype = float) - np.array(a_pos[0], dtype = float)
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