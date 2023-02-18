from IPython.display import clear_output
import numpy as np
import fractions as f
import math
from scipy.spatial import Voronoi, voronoi_plot_2d, distance

def update_progress(progress):
  
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
        
    block = int(round(bar_length * progress))

    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "■" * block + "○" * (bar_length - block), progress * 100)
    print(text)

def cart2xys(pos, cell):
    """
    Cartesian (absolute) position in angstrom to fractional position (scaled position in lattice).
    """
    pos = np.asarray(pos)
    bg = np.linalg.inv(cell.lattice)
    xyzs = np.tensordot(bg, pos.T, axes=([-1], 0)).T
    return xyzs

def xys2cart(pos, cell):
    """
    Fractional position (scaled position in lattice) to cartesian (absolute) position in angstrom.
    """
    pos = np.asarray(pos)
    xyzr = np.tensordot(cell.lattice, pos.T, axes=([-1], 0)).T
    return xyzr

def translate_pretty(fractional, pbc):
    """Translates atoms such that fractional positions are minimized."""

    for i in range(3):
        if not pbc[i]:
            continue

        indices = np.argsort(fractional[:, i])
        sp = fractional[indices, i]

        widths = (np.roll(sp, 1) - sp) % 1.0
        fractional[:, i] -= sp[np.argmin(widths)]
        fractional[:, i] %= 1.0
    return fractional
    
# def to_ase_atoms()

def to_sublists(lst, length=2):
    """
    Split a list to sublists with specified length: e.g. a = [a,b,c,d,e]
    to_sublists(a) => [[a,b], [b,c], [c,d], [d,e]] 
    """
    return [lst[i:i+length] for i in range(len(lst)+1-length)]

def retrieve_symbol(string):
    """function to remove numbers in a string, so that the atom dict keys can be converted to chemical symbols"""
    return ''.join([i for i in string if not i.isdigit()])
  
def get_gcd(ints):
    """
    Calculate the maximal common divisor of a list of integers
    """
    gcd = math.gcd(ints[0],ints[1])
    for i in range(2,len(ints)):
        gcd = math.gcd(gcd,ints[i])
    return gcd


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
        # if self.name == "m" or self.name == "miller":
        #     self.array = self.array.astype(int)
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
        # print(self.array)
        return self
        
    def to_cart(self, normalise = True):
        if self.name == "c" or self.name == "cartesian":
            print("Already cartesian!")
        else:
            self.array = np.dot(self.array.T, self.cell)
            if normalise:
                self.array = self.array / np.linalg.norm(self.array)
            self.name = "cartesian"
        # print(self.array)
        return self


class h_channel(Voronoi):
    def __init__(self, points, furthest_site=False, incremental=False, qhull_options=None):
        super().__init__(points, furthest_site, incremental, qhull_options)

    def sort_points(self, points):
        indices = [self.point_region[distance.cdist([point], self.points).argmin()] for point in points]
        return indices

def calc_freq(regions, timestep = None):
    """
    Calculates the frequency of switching channels, unit times per ps
    """
    count = 0
    current_region = regions[0]
    freq = None
    
    if timestep is None:
        print("No timestep specified, defauting to 5 fs!")
        timestep = 5
        
    for i, region in enumerate(regions):
        if region != current_region:
            count += 1
            current_region = region
        else:
            pass
        
        freq = count * 5/ ((i + 1) * 1000)
        
    return freq

def get_basis(h_path, cell, miller = True):
    """Identifies two vectors that are perpendicular to the h_path,
    if h_path is in Cartesian coordinates rather than Miller indices, set miller = False"""
    # Create the vector object
    if miller:
        h_path = vector(h_path, cell = cell, name = "m")
        h_path = h_path.to_cart(normalise = True)
    else:
        h_path = vector(h_path, cell = cell, name = "c")
        h_path = h_path.to_cart(normalise = True)
    
    # Specify one vector that is orthogonal with the h_path
    indices = [0, 1, 2]
    idx = np.where(h_path.array != 0)[0][0]
    indices.remove(idx)
    basis_x = np.ones(3)
    basis_x[idx] = -(h_path.array[indices[0]] + h_path.array[indices[1]])/h_path.array[idx]
    basis_x = basis_x/np.linalg.norm(basis_x)
    
    # Calculate the other vector as the cross product
    basis_y = np.cross(h_path.array, basis_x)
    return basis_x, basis_y

def trans_coord(points, cell, W, miller = True):
    """Change coordinate system of points from [[1,0,0],[0,1,0],[0,0,1]] to W"""
    # We need to get cartesion coordinates if the points are in Miller indices
    trans_pos = np.zeros(points.shape)
    if miller:
        for i, point in enumerate(points):
            point = vector(point, cell).to_cart(normalise = False).array
            trans_pos[i] = np.linalg.inv(W).dot(point)
    return trans_pos
