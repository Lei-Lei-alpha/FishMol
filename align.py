import numpy as np
import pandas as pd

def TFA_align(frames, at_dicts, t0 = 0, filename = None):
    """
    Filter the trajctory and save specified atoms only.
    contents: the trajectory file.
    """
    columns=["t"]
    for i in range(1,len(at_dicts)+1):
        columns += [f"TFA{i}_x", f"TFA{i}_y", f"TFA{i}_z"]
    align_df = pd.DataFrame(columns=columns)
    
    t_end = t0 + 5*(len(frames)-1) # 5 fs is the interval of traj
    align_df["t"] = np.linspace(t0, t_end, num = len(frames))
    
    for i, frame in enumerate(frames):
        for j, at_dict in enumerate(at_dicts):
            indices = [at_dict["C1"], at_dict["C2"]]
            positions = [frame[i,0]["position"] for i in indices]
            alignment = [pos for pos in positions[1] - positions[0]]
            align_df.loc[i, [f"TFA{j+1}_x", f"TFA{j+1}_y", f"TFA{j+1}_z"]] = alignment
        
    if filename != None:
        align_df.to_excel(filename)
    
    return align_df

def theta(vec1, vec2):
    """
    Calculate the angle between two vectors.
    """
    return np.arccos(vec1.dot(vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))*180/np.pi

def reorient_theta(data):
    """
    Calculates the alignment of C-C bond to current alighnment. The alighnment changes when the C-C bond direction changed for more than 80 degrees.
    """
    orientation = data[0]
    alignment = data[0]
    output = []
    reorientations = []
    flips = []
    for i in range(len(data)):
        angle = theta(np.array(data[i]), np.array(orientation))
        tot_rot = theta(np.array(data[i]), np.array(alignment))
        rotate = angle >= 80
        flip = tot_rot >= 180
        if rotate:
            if flip:
                alignment = data[i]
                flips.append([i*5/1000, angle])
                print(f"Rotated & Flipped!\n Alignment changed to {alignment},\n current frame {i}, \n orientation changed to {orientation},\n current frame {i}, \n current time {i*5/1000} ps, \n rotation from previous orientation {angle}, \n rotation from previous alignment {tot_rot}")
      
            orientation = data[i]
            reorientations.append([i*5/1000, angle])
            print(f"Rotated!\n Orientation changed to {orientation},\n current frame {i}, \n current time {i*5/1000} ps, \n rotation from current alignment {tot_rot}, \n rotation from previous orientation {angle}")
            
        elif flip:
            alignment = data[i]
            flips.append([i*5/1000, angle])
            print(f"Flipped!\n Alignment changed to {alignment},\n current frame {i}, \n current frame {i}, \n current time {i*5/1000} ps, \n rotation from previous orientation {angle}, \n rotation from previous alignment {tot_rot}")
                
        output.append([angle, list(orientation), list(alignment)])
    return reorientations, flips, output