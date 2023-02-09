import numpy as np
import pandas as pd
from fishmol import c1_at_dict, msd

# Define functions to convert vectors between miller indices and cartesian coordinates
import fractions as f
import math
def get_gcd(ints):
    """
    Calculate the maximal common divisor of a list of integers
    """
    gcd = math.gcd(ints[0],ints[1])
    for i in range(2,len(ints)):
        gcd = math.gcd(gcd,ints[i])
    return gcd

class vector:
    """
    Vecotor object that can convert between miller indices and cartesian coordination (normalised).
    """
    def __init__(self, array, cell, name = "miller"):
        if len(array) != 3:
            raise Exception("The input array must has a length of exactly 3.")
        else:
            self.array = np.asarray(array)
        if name not in ["m", "miller", "c", "cartesian"]:
            raise Exception("Unrecognised name! Please use 'm' or 'miller' if the input array is a miller index, use 'c' or 'cartesian' if it is a coordinate.")
        else:
            self.name = name
        if self.name == "m" or self.name == "miller":
            self.array = self.array.astype(int)
        self.cell = np.asarray(cell)
        
    def to_miller(self):
        if self.name == "m" or self.name == "miller":
            print("Already miller!")
        else:
            # obtain coord refer to the lattice vector
            if len(self.cell) == 3:
                pass
                
            elif len(self.cell) == 6:
                self.array = np.dot(self.array, self.cell.T)
                
            h = f.Fraction(self.array[0]).denominator
            k = f.Fraction(self.array[1]).denominator
            l = f.Fraction(self.array[2]).denominator
            
            self.array = self.array * h * k * l

            self.array = self.array.astype(int) // get_gcd(self.array.astype(int))
            self.name = "miller"
        print(self.array)
        
    def to_cart(self, normalise = True):
        if self.name == "c" or self.name == "cartesian":
            print("Already cartesian!")
        else:
            self.array = np.dot(self.array.T, self.cell)
            if normalise:
                self.array = self.array / np.linalg.norm(self.array)
            self.name = "cartesian"
        print(self.array)
        

phi = np.linspace(-np.pi, np.pi, num = 20)
theta = np.linspace(-np.pi/2, np.pi, num = 20)

x = np.sin(phi) * np.sin(theta)
y = np.sin(phi) * np.cos(theta)
z = np.cos(phi)

X, Y, Z = np.meshgrid(x, y, z)

vecs = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis = 1)

ave_D = np.zeros(len(vecs))

t0 = 0
t_end = t0 + 5*(len(water_com)-1) # 5 fs is the interval of traj
t = np.linspace(t0, t_end, num = len(water_com))

t = t[1001:]
t -= 5000
    
water_com = pd.read_excel("cage1/water_com_400K.xlsx", header=0, index_col=0, engine = "openpyxl")

cell= [
    [21.2944000000,        0.0000000000,        0.0000000000],
    [-4.6030371123,       20.7909480472,        0.0000000000],
    [-0.9719093466,       -1.2106211379,       15.1054299403]
]


for i, vec in enumerate(vecs):
    msd_df = np.zeros((len(water_com) - 1, 8))
    for j in range(8):
        temp = np.array(water_com.iloc[:, 3*j:3*j+3])
        temp = traj_proj(temp, vec)
        msds = msd.msd_1d(temp)[1:]
        msd_df[:, j] = msds

    ave_msd = msd_df.mean(axis=1)
    ave_msd = ave_msd[1000:]

    linear_model=np.polyfit(t/1000, ave_msd, 1)
    ave_D[i] = linear_model[0]/20000
    
results = pd.DataFrame()
results[["x", "y", "z"]] = vecs
results["mean_D"] = ave_D
results.to_excel("cage1/400K_anisotropy_all_direc.xlsx")