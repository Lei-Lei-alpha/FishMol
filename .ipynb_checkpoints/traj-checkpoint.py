#!/usr/bin/env python
from ase import Atoms
import pandas as pd
import numpy as np
import mmap
from .data import elements

def read(traj, prop = None):
    """
    read the xyz trajctory file, and remove the headers: line1: number of atoms in the system, line2: Properties of the system, pbc, simu step, time, a,b,c etc.
    prop: input the start of second header if it is not "Properties"
    """
    with open(traj) as f:
        mm = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
        frames = [line.rstrip().decode('utf8') for line in iter(mm.readline, b"")]
        mm.close()
    f.close()
    header = frames[0]
    natoms = int(header)
    nframes = frames.count(frames[0])
    if prop == None:
        prop = frames[1].split("=")
    dt = np.dtype([('symbol', np.unicode_, 16), ('position', np.float64, (3,))]) # numpy datatype object, one element is string , the other is 3 floats
    # Store each line as a numpy array
    frames = [np.array([(line.split(" ")[0],
                         line.split(" ")[-3:])],
                       dtype=dt) 
              for line in frames if not (line.startswith(header) or line.startswith(prop[0]))]  
    # split the trajectory into frames
    frames = [np.array(frames[x: x + natoms]) for x in range(0, len(frames), natoms)]
    return natoms, nframes, frames

def retrieve_symbol(string):
    """function to remove numbers in a string, so that the atom dict keys can be converted to chemical symbols"""
    return ''.join([i for i in string if not i.isdigit()])

def get_com(frame, at_dict):
    """
    Calculate the center of mass 
    """
    indices = list(at_dict.values())
    symbols = [frame[i,0]["symbol"] for i in indices]
    positions = [frame[i,0]["position"] for i in indices]
    masses = np.array([elements[symb] for symb in symbols])
    com = np.dot(masses, positions)/masses.sum()
    return com

def calc_com(frames, at_dicts, filename = "water_com.xlsx"):
    """
    Calculate the center of mass of water molecules and write the results into a excel file.
    """
    # make pandas.DataFrame to store data
    columns=[]
    for i in range(1,len(at_dicts)+1):
        columns += [f"water{i}_x", f"water{i}_y", f"water{i}_z"]
    df_water_com = pd.DataFrame(columns=columns)
    
    for frame in frames:
        # Make an empty list to store centre of masses
        com = []
        for at_dict in at_dicts:
            # Calculate com
            com += [x for x in get_com(frame, at_dict)]
        # Add data to dateframe
        df_water_com.loc[len(df_water_com.water1_x)]=com
    
    # Write data to excel file
    df_water_com.to_excel(filename)
    print("Done! CoM data wrote to " + filename)
    return df_water_com

def filter_traj(frames, at_dicts, filename = None):
    """
    Filter the trajctory and save specified atoms only.
    contents: the trajectory file.
    """
    for frame in frames:
        system = Atoms('', positions=[])
        for at_dict in at_dicts:
            indices = list(at_dict.values())
            symbols = [frame[i,0]["symbol"] for i in indices]
            positions = [frame[i,0]["position"] for i in indices]

            atom = Atoms(symbols,
                        positions = positions,
                       )
            system += atom
        system.write(filename, append = True)
