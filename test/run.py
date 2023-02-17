import subprocess
cages = [1,2]
temps=[400,500,600]
for cage in cages:
    for temp in temps:
        print(f"AIMD results Cage-{cage} {temp}K:")
        subprocess.call(f"python3 md_analysis.py -cage {cage} -temp {temp}", shell=True)
        # subprocess.call(f"python3 hb_analysis.py -cage {cage} -temp {temp}", shell=True)