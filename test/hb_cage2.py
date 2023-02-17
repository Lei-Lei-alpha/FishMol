print("Starting ... ...")

import mmap
import numpy as np

def read_traj(traj, prop = None, select = ":"):
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
    
    if prop == None:
        prop = contents[1].split("=")
        
    dt = np.dtype([('symbol', np.unicode_, 2), ('position', np.float64, (3,))])
    
    if select == ":" or select == "all":
            nframes = contents.count(contents[0])
        
    if isinstance(select, slice):
        start = select.start
        stop = select.stop
        step = select.step
        if start is None:
            start = 0
        if stop is None:
            stop = contents.count(contents[0])
        if step is None:
            step = 1
        contents = contents[start*(natoms+2):stop*(natoms+2):step]
        nframes = contents.count(contents[0])
    
    contents = [np.asarray((line.split()[0], line.split()[-3:]), dtype=dt) 
              for line in contents if not (line.startswith(header) or line.startswith(prop[0]))]
    
    contents = [np.asarray(contents[x: x + natoms]) for x in range(0, len(contents), natoms)]
    
    return natoms, nframes, contents

print("Reading trajectory file ... ...")
natoms, nframes, frames = read_traj("cage2-500K.xyz", select = slice(1000,11001,1))

from md_analysis import hbond
from md_analysis import c2_at_dict

def pop_C(TFAs):
    for TFA in TFAs:
        [TFA.pop(key) for key in list(TFA.keys()) if key == "C1" or key == "C2"] # keep F and O in TFAs only for H-bond analysis
    return TFAs

waters = c2_at_dict.waters
TFA_fracs = c2_at_dict.TFAs
TFAs = c2_at_dict.TFAs
phenols = c2_at_dict.phenols
amines = c2_at_dict.amines
TFA_F_Os = pop_C(TFAs)
atoms_dict = waters + TFA_fracs + phenols + amines + TFA_F_Os

donors, acceptors = hbond.find_d_a(atoms_dict = atoms_dict)

d_a_pairs, d_h_a_pairs = hbond.pair_d_a(donors = donors, acceptors = acceptors)

cell= [
    [16.5370000000, 0.0000000000, 0.0000000000],
    [1.0363626153, 20.3876763887, 0.0000000000],
    [4.1539399394, 3.1417959791, 20.6513608512]
]

import itertools
d_a_pairs = [list(itertools.chain.from_iterable([d_a_pair[0].values(), d_a_pair[1].values()])) for d_a_pair in d_a_pairs]
d_h_a_pairs = [list(itertools.chain.from_iterable([d_h_a_pair[0].values(), d_h_a_pair[1].values()])) for d_h_a_pair in d_h_a_pairs]

from md_analysis.data import vdW_R
def cart2xys(pos, cell):
    """
    Cartesian (absolute) position in angstrom to fractional position (scaled position in lattice).
    """
    pos = np.asarray(pos)
    bg = np.linalg.inv(cell)
    xyzs = np.tensordot(bg, pos.T, axes=([-1], 0)).T
    return xyzs

def xys2cart(pos, cell):
    """
    Fractional position (scaled position in lattice) to cartesian (absolute) position in angstrom.
    """
    pos = np.asarray(pos)
    xyzr = np.tensordot(cell, pos.T, axes=([-1], 0)).T
    return xyzr

def theta(vec1, vec2):
    """
    Calculate the angle between two vectors.
    """
    return np.arccos(np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))*180/np.pi

def distance(atom_pair, frame, cell = None, mic = False):
    """
    Calculates the distance of a atom pair.
    """
    positions = frame[atom_pair]["position"]

    if mic:
        a2b = cart2xys(positions[1], cell = cell) - cart2xys(positions[0], cell = cell)
        for i in range(3):
            a2b[i] = a2b[i] - round(a2b[i])
        a2b = xys2cart(a2b, cell = cell)
    else:
        a2b = positions[1] - positions[0]
    return np.sqrt(np.dot(a2b, a2b))

def res_h(d_a_pairs, d_h_a_pairs, frames, cell = None):
    """
    Calculates the H-bond information in each frame of the trajectory.
    """
    results = []
    for i, frame in enumerate(frames):
        hbonds = [] #list to store all hbonds in current frame
        for d_a_pair, d_h_a_pair in zip(d_a_pairs, d_h_a_pairs):
            # Sum of van der Waals radii
            symbols = frame[d_a_pair]["symbol"]
            vdW_sum = vdW_R[symbols[0]] + vdW_R[symbols[1]]
            
            # not a H-bond of D-A distance is greater than their vdW radii times 1.05, 1.05 to take bond length change during MD simulation.
            r_d_a = distance(d_a_pair, frame, cell = cell, mic = True)  # calculate the D-A distance
            if r_d_a < 1.02 * vdW_sum:
                # calculate the D-H⋅⋅⋅A angle
                d_h_pos = frame[d_h_a_pair[:2]]["position"] # the positions of donor and hydrogen
                a_pos = frame[d_h_a_pair[2]]["position"] # the positions of acceptor
                d_h_vec = d_h_pos[1] - d_h_pos[0]
                a_h_vec = d_h_pos[1] - a_pos
                d_h_a_ang = theta(d_h_vec, a_h_vec) # angle
                d_h = distance(d_h_a_pair[:2], frame, cell = cell, mic = True) # calculate the D-H length
            
                # the D-H⋅⋅⋅A angle criteria used: the D-H⋅⋅⋅A angle is close to a right angle refer to the D-H⋅⋅⋅A angle - R(D⋅⋅⋅A) plot
                # an angle range is included considering the oscillation of bond lenghth and anlgle
                if d_h_a_ang >= (np.rad2deg(np.arctan2(r_d_a, d_h)) + 180)*3/8:
                # if d_h_a_ang >= 90:
                    # Store current H-bond
                    hbonds.append(
                          {
                              "donor": [frame[d_h_a_pair[:2]]["symbol"].tolist(), d_h_a_pair[:2]],
                              "acceptor": [frame[d_h_a_pair[2]]["symbol"].tolist(), d_h_a_pair[2]],
                              "R(D-A)": r_d_a,
                              "DHA_ang": d_h_a_ang,
                              "r(D-H)": d_h,
                          }
                    )
                
        results.append({f"frame": i, "n_hbonds": len(hbonds), "hbonds": hbonds})
    return results

print("Resolving H-bonds ... ...")

hbonds = res_h(d_a_pairs = d_a_pairs, d_h_a_pairs = d_h_a_pairs, frames = frames, cell = cell)
import json
with open("hbonds_mic_c2_500K.json", "w") as json_file:
    json.dump(hbonds, json_file)
json_file.close()

print("Done!")