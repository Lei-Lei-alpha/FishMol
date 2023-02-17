"""
Run the analysis by passing arguments, e.g. python3 hb_analysis.py -cage 1 -tem 400
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

from md_analysis import hbond
from md_analysis import c1_at_dict
from md_analysis import c2_at_dict
import argparse
import os
import json

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

def pop_C(TFAs):
    for TFA in TFAs:
        [TFA.pop(key) for key in list(TFA.keys()) if key == "C1" or key == "C2"] # keep F and O in TFAs only for H-bond analysis
    return TFAs

# Read Trajectory file
print(f"Reading traj file {cage}-{temperature}.xyz ...")
_, _, frames = hbond.read_traj(f"{cage}-{temperature}.xyz")

# Atom_dicts here
if cage == "cage1": # Cage 1 atom_dicts
    waters = c1_at_dict.waters
    TFAs = c1_at_dict.TFAs
    phenols = c1_at_dict.phenols
    amines = c1_at_dict.amines
    TFA_F_Os = pop_C(TFAs)
    atoms_dict = waters + phenols + amines + TFA_F_Os

else:
    waters = c2_at_dict.waters
    TFA_fracs = c2_at_dict.TFAs
    TFAs = c2_at_dict.TFAs
    phenols = c2_at_dict.phenols
    amines = c2_at_dict.amines
    TFA_F_Os = pop_C(TFAs)
    atoms_dict = waters + TFA_fracs + phenols + amines + TFA_F_Os

donors, acceptors = hbond.find_d_a(atoms_dict = atoms_dict)
d_a_pairs, d_h_a_pairs = hbond.pair_d_a(donors = donors, acceptors = acceptors)

print("Analysing all H-bonds...")

hbonds = hbond.res_h(d_a_pairs=d_a_pairs, d_h_a_pairs=d_h_a_pairs, frames = frames)

with open(outdir + f"hbonds_{temperature}.json", "w") as json_file:
    json.dump(hbonds, json_file)
json_file.close()

# Analysis of H-bonds involved with water molecules only
water_d, water_a = hbond.find_d_a(atoms_dict = waters)

print("Done! Analysing H-bonds formed with waters only...")
hbonds = hbond.res_h(d_a_pairs=d_a_pairs, d_h_a_pairs=d_h_a_pairs, frames = frames, water_d = water_d, water_a = water_a)
with open(outdir + f"hbonds_water_{temperature}.json", "w") as json_file:
    json.dump(hbonds, json_file)
json_file.close()

print(f"All done! Results saved to {outdir}")