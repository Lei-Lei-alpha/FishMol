import pandas as pd
import numpy as np

def get_com(frame, at_dict):
    """
    Calculate the center of mass 
    """
    indices = list(at_dict.values())
    symbols = [frame[i,0]["symbol"] for i in indices]
    positions = [frame[i,0]["position"] for i in indices]
    masses = np.array([data.elements[symb] for symb in symbols])
    com = np.dot(masses, positions)/masses.sum()
    return com

def calc_com(frame, at_dicts, filename = "water_com.xlsx"):
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
        for at_dict in enumerate(at_dicts):
            # Calculate com
            com += [x for x in get_com(frame, at_dict)]
        # Add data to dateframe
        df_water_com.loc[len(df_water_com.water1_x)]=com
    
    # Write data to excel file
    df_water_com.to_excel(filename)
    print("Done! CoM data wrote to " + filename)
    return None