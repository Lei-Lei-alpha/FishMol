"""
Run the analysis by passing arguments, e.g. python3 md_analysis.py -cage 1 -temp 400
"""

print("#  ○                     ○                       ▄▄▄▀\n\
#                                            ▄▄███▀\n\
#               ○                         ▄▄█████▀          FishMol version 0.0\n\
#     ○                                ▄▄████████▄\n\
#                ○                  ▄▄████████████▄                  ○\n\
#                            ▄▄▄████████████████████▄\n\
#   ○                  ▄▄▄██████ ██████████████████████▄▄\n\
#          ○       ▄▄████████████ █████████████████████████▄▄                      ○   ▄▄▀\n\
#               ▄█████████████████ ████████████████████████████▄▄                   ▄███▀\n\
#        ○   ▄██████████▀▀▀████████ ██████████████████████████████▄▄             ▄█████▀\n\
#           ▄████████▀   □  ▀███████ █████████████████████████████████▄▄      ▄▄██████▀\n\
#           ■▄███████▄      ▄██████ ██████████████████████████████████████▄▄▄████████\n\
#            ▀█████████▄▄▄████████ █████████████████████████████████████████████████\n\
#              ▀█████████████████ ██████████████████████████████████████████████████▄\n\
#                ▀█████████████ █████████████████████████████████████▀▀▀▀▀████████████▄\n\
# ○                ▀█████████ █████████████████████████████████▀▀           ▀▀████████▀\n\
#                     ▀▀███████████████████████████████▀▀▀▀██▀                 ▀▀████▀\n\
#                          ▀▀▀█████████████████▀▀▀▀         ▀▄     ○              ▀█▀\n\
#                                ▀▀███████▀▀\n\
#                 ○                 ▀████\n\
#                                      ▀▀▄")
print("Starting...")

import numpy as np
import pandas as pd
from md_analysis import traj
from md_analysis import msd
from md_analysis import plot
from md_analysis import align
from md_analysis import c1_at_dict
from md_analysis import c2_at_dict
import json
import os
import argparse

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

# Read Trajectory file
print(f"Reading traj file {cage}-{temperature}.xyz ...")
natoms, nframes, frames = traj.read(f"{cage}-{temperature}.xyz")

# Atom_dicts here
if cage == "cage1": # Cage 1 atom_dicts
    waters = c1_at_dict.waters
    TFAs = c1_at_dict.TFAs
    # Anisotropy of diffusion
    vecs ={ 
        "x": [1, 0, 0], 
        "y": [0, 1, 0], 
        "z": [0, 0, 1], 
        "path": [1, 1, 0],
    }

else: # Cage 2 atom_dicts
    waters = c2_at_dict.waters
    TFAs = c2_at_dict.TFAs
    # Anisotropy of diffusion
    vecs ={ 
        "x": [1, 0, 0], 
        "y": [0, 1, 0], 
        "z": [0, 0, 1], 
        "path": [1, 0, 1],
    }

# Calculate CoM of water molecules

# Params
header = 2
formula = "OH2"

print("Calculating the cetre of mass for water molecules...")

water_com = traj.calc_com(frames = frames, at_dicts = waters, filename = outdir + f"water_com_{temperature}.xlsx")

# Calculate the MSD and D
print("Calculating the MSDs and diffusion coefficients...")

t0 = 0
t_end = t0 + 5*(nframes-1) # 5 fs is the interval of traj
t = np.linspace(t0, t_end, num = nframes)

columns = ["t"]

for i in range (1, len(waters) + 1):
    columns += [f"water{i}_MSD", f"water{i}_D"]

columns += ["Mean_MSD", "MSD_error", "Mean_D", "D_error"]

msd_d_df = pd.DataFrame(columns = columns)

msd_d_df["t"] = t[1:]

for i in range(len(waters)):
    df = np.array(water_com.iloc[:, 3*i:3*i+3])
    msds = msd.msd_fft(df)
    msd_d_df[columns[i*2 + 1]] = msds[1:]
    msd_d_df[columns[i*2 + 2]] = msds[1:] * 1E-16 / (6*t[1:]*1E-15)

msd_df = msd_d_df.iloc[:, 1:2*len(waters):2]
d_df = msd_d_df.iloc[:, 2:2*len(waters)+1:2]

msd_d_df["Mean_MSD"] = msd_df.mean(axis=1)
msd_d_df["MSD_error"] = msd_df.std(axis=1) / len(waters)**0.5

msd_d_df["Mean_D"] = d_df.mean(axis=1)
msd_d_df["D_error"] = d_df.std(axis=1) / len(waters)**0.5

# Save data
msd_d_df.to_excel(outdir + f"MSD_D_H2O_{temperature}.xlsx")

# Plot and save figure
# Average
plot.ave_msd_d(msd_d_df, filename = outdir + f"ave_msd_d_{temperature}.jpg")
# Molecule compare
plot.msd_d_comp(msd_d_df, waters, filename = outdir + f"msd_d_comp_{temperature}.jpg")

print("Calculating the anistropy of MSDs and diffusion coefficients...")

# ani_df = pd.DataFrame()
# ani_df["t"] = t[1:]

# for i in range (len(waters)):
#     for name, vec in vecs.items():
#         temp_com = water_com.iloc[:,3*i:3*i+3]
#         temp_com = msd.traj_proj(temp_com, vec)
#         temp_msd = msd.msd_1d(temp_com)[1:]
#         ani_df[f"water{i+1}_{name}_MSD"] = temp_msd
#         ani_df[f"water{i+1}_{name}_D"] = temp_msd* 1E-16 / (2*t[1:]*1E-15)

# ani_df.to_excel(outdir + f'anisotropy{temperature}.xlsx')

# direction = list(vecs.keys())

# for i in range(1,len(waters)+1):
#     index = i
#     filename = outdir + f"water{i}.jpg"
#     plot.plot_comp(df=ani_df, index = index, direction = direction, filename = filename)

# TFA alignment

print("Resolving the alignment of TFA anions...")

align_df = align.TFA_align(frames, TFAs, filename = outdir + f"TFA_CC_align_{temperature}.xlsx")

results = {}
for i in range(len(TFAs)):
    print(f"TFA{i+1}:")
    reorientations, flips, output = align.reorient_theta(np.array(align_df.iloc[:,3*i+1:3*i+4]))
    plot.quick_plot(align_df.iloc[:,0],
    np.array(output, dtype=object)[:,0],
    reorientations = reorientations,
    flips = flips,
    filename = outdir + f"TFAs{i+1}_{temperature}.jpg")
    results[f"TFA{i+1}"]={"reorientations": reorientations,
                          "flips": flips,
                          "orientations": output,}
with open(outdir + f"TFAs_align_{temperature}.json", "w") as json_file:
    json.dump(results, json_file)
json_file.close()

print(f"All done! Results saved to {outdir}")