#! /usr/bin/python3

import os
import subprocess
import argparse
import pandas as pd
from fishmol import trj
from fishmol import atoms
from fishmol import utils
from fishmol import funcs
from fishmol.cages import c1_at_dict, c2_at_dict

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

plt.rcParams.update({
    "font.size": 13,
    "font.family": "sans-serif",
    # "font.sans-serif": "Arial",
    # "font.weight": "heavy",
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "figure.figsize": (4.2,3.6),
    "figure.subplot.left": 0.21,
    "figure.subplot.right": 0.96,
    "figure.subplot.bottom": 0.18,
    "figure.subplot.top": 0.93,
    "legend.frameon": False,
})

from colour import Color
from matplotlib.colors import LinearSegmentedColormap
ramp_colors = ["#ffffff", "#9ecae1", "#2166ac", "#1a9850", "#ffff33", "#b2182b", "#67000d"]
color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )

# Check if fishmol is installed
subprocess.run("fishmol")

# Construct an argument parser
all_args = argparse.ArgumentParser()

# Add arguments to the parser
all_args.add_argument("-cage", "--Value1", required=True,
   help="first Value")
all_args.add_argument("-temp", "--Value2", required=True,
   help="second Value")
args = vars(all_args.parse_args())

cage = "cage" + args['Value1']

temperature = args['Value2'] + "K"

outdir = f"results/{cage}/nvt-{temperature}/"

if not os.path.exists(outdir):
    print(f"{outdir} does not exist, creating...")
    os.makedirs(outdir)


# Atom_dicts here
if cage == "cage1": # Cage 1 atom_dicts
    waters = c1_at_dict.waters
    TFAs = c1_at_dict.TFAs
    # Anisotropy of diffusion
    cell = [
        [21.2944000000,        0.0000000000,        0.0000000000],
        [-4.6030371123,       20.7909480472,        0.0000000000],
        [-0.9719093466,       -1.2106211379,       15.1054299403]
    ]

else: # Cage 2 atom_dicts
    waters = c2_at_dict.waters
    TFAs = c2_at_dict.TFAs
    # Anisotropy of diffusion
    cell = [
        [16.5370000000, 0.0000000000, 0.0000000000],
        [1.0363626153, 20.3876763887, 0.0000000000],
        [4.1539399394, 3.1417959791, 20.6513608512]
    ]

# Load data
print(f"Reading traj file {cage}-{temperature}.xyz ...")
traj = trj.Trajectory(timestep = 5, data = f"{cage}-{temperature}.xyz", index = slice(1000, None, 1), cell = cell)

# Water VRD

print("Calculating the VRD of water molecules ...")
water_vrd = pd.DataFrame()

for i, water in enumerate(waters):
    spec = [*water.values()]
    vrd = funcs.VRD(traj = traj, spec = [[spec[0],], spec[1:]], num = 500, sampling = 10, skip = 2)
    results = vrd.calculate(plot = False)
    water_vrd[f"t{i+1}"] = pd.Series(results.t)
    water_vrd[f"w{i+1}_u"] = pd.Series(results.c_t_mean[:,0])
    water_vrd[f"w{i+1}_u_e"] = pd.Series(results.c_t_error[:,0])
    water_vrd[f"w{i+1}_v"] = pd.Series(results.c_t_mean[:,1])
    water_vrd[f"w{i+1}_v_e"] = pd.Series(results.c_t_error[:,1])

water_vrd.to_excel(outdir + f"water_vrd_{temperature}.xlsx")

# TFA VRD

print("Calculating the VRD of TFA molecules ...")
tfa_vrd = pd.DataFrame()

for i, TFA in enumerate(TFAs):
    spec = [*TFA.values()]
    vrd = funcs.VRD(traj = traj, spec = [[spec[-2],], [spec[-1],]], num = 5000, sampling = 10, skip = 5)
    results = vrd.calculate(plot = False)
    tfa_vrd[f"t{i+1}"] = pd.Series(results.t)
    tfa_vrd[f"w{i+1}_u"] = pd.Series(results.c_t_mean[:,0])
    tfa_vrd[f"w{i+1}_u_e"] = pd.Series(results.c_t_error[:,0])

tfa_vrd.to_excel(outdir + f"TFA_vrd_{temperature}.xlsx")