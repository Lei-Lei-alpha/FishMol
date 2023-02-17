"""
Run the analysis by passing arguments, e.g. python3 md_analysis.py -cage 1 -temp 400
"""

print(
    "################################################################################\n\
#####  #####     ###  ##########===============###########  #####     ###  #####\n\
#####  #####  ######  #########|| MD Analysis ||##########  #####  ######  #####\n\
#####  #####     ###  #########|| version 0.0 ||##########  #####     ###  #####\n\
#####  #####  ######  #########|| MIT licence ||##########  #####  ######  #####\n\
#####     ##     ###  ##########===============###########     ##     ###  #####\n\
################################################################################\n"
    )

print("Starting...")

import numpy as np
from md_analysis import traj
from md_analysis import plot
from md_analysis import align
from md_analysis import c1_at_dict
from md_analysis import c2_at_dict
import pandas as pd
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


if cage == "cage1": # Cage 1 atom_dicts
    TFAs = c1_at_dict.TFAs

else: # Cage 2 atom_dicts
    TFAs = c2_at_dict.TFAs

# TFA alignment

print("Resolving the alignment of TFA anions...")

align_df = pd.read_excel(outdir + f"TFA_CC_align_{temperature}.xlsx", header=0, index_col=0, engine ="openpyxl")

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