import re
import json
import warnings
import itertools
import numpy as np
from fishmol import trj, data
from recordclass import  dataobject, astuple, asdict
from iteration_utilities import deepflatten
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from colour import Color
from matplotlib.colors import LinearSegmentedColormap
from fishmol.utils import to_sublists, update_progress, mic_dist

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

class dimers(object):
    """The dimers in each frame of a trajectory
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Arguments
    at_g1 : list, the indices of the donors of H-bonds or one of the molecules that forms the dimer. For the case of H-bond donors, pass the indices as a list in the following format[D, [H1, H2, ...]] (eg:  for two water molecules whose indices are O: 0, H: 1, H: 2; and O: 8, H: 8, H: 10, the at_g1 should be either [[0, [1, 2]], [8, [9, 10]]].
    at_g2 : list or tuple, the indices of the acceptors or the ther molcecule the forms the dimer.
    criteria : list of geometry criteria for the specified atoms to be regarded as pairs. Accepted format for distance:
        1. number + vdw_sum or number + coval_sum. eg: 1.05vdw_sum or 1.1coval_sum, which keep all pairs whose distances are smaller than their sum of van der Waals radii times 1.05 or their sum of covalency radii times 1.1. Plese set com = False when using this type of distance criteria.
        2. Number specifying the distance of two atoms (or the distance of the centre of masses of two atom groups by setting com = True) in angstroms.
        For angles, pass the angle of donor-H⋅⋅⋅acceptor (∠DHA). It also possible to use a flexible angle criteria to take the dynamics into consideration. To do this, pass a string containing a number of specifying the desired angle in degree and the letter "f", eg: "120f".
    com : wether to use the centre of mass of the atom groups in the distance calculation. Set com = True if the atoms groups are two lists of molecules.
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Methods:
    make_pairs: prepare the unique pairs from given atom groups for further filter based on given criteria. Returns None, the pairs are stored in results.pairs.
    pair_filter: filter the pairs and find dimers in each frame of the trajectory according to the specified criteria. Returns None, the dimers in each frame of the trajectory are stored in results.dimers. If a filename is specified, the resolved dimer data will be writen in the give file name in .json format.
    """
    def __init__(self, traj, at_g1, at_g2, criteria = [None, None], com = False):
        self.traj = traj
        self.g1 = at_g1
        self.g2 = at_g2
        criteria_dc = make_dataclass("criteria", "distance angle")
        self.criteria = criteria_dc(*criteria)
        if self.criteria.distance is None:
            warnings.warn("Distance criteria is mandatory! The dimer will be resolved by default distance settings! To ensure the accuracy of dimer resolvation, specify the geometry criteria by passing a distance (and an angle if applicable)")
            self.criteria.distance = "1.05vdw_sum"
        self.com = com
        results = make_dataclass("Results", "pairs dimers")
        self.results = results
                  
    def make_pairs(self):
        """
        Make d_a_pairs and d_h_a_pairs from the donors and acceptors dict.
        """
        # com = self.com
        # When calculates the distance of centre of masses of two molecules groups
        if self.com:
            self.results.pairs = list(map(list, [x for x in itertools.product(self.g1, self.g2) if x[0] != x[1]]))
            # Remove any duplicates if at_g1 and at_g2 are the same
            for pair1 in self.results.pairs:
                for pair2 in self.results.pairs:
                    if pair1[0] == pair2[1] and pair1[1] == pair2[0]:
                        self.results.pairs.remove(pair2)
            self.results.pairs = np.asarray(self.results.pairs)
        # Otherwise
        else:
            # simple case: make pairs from the the listed atoms in at_g1 and at_g2
            if self.criteria.angle is None:
                if all([isinstance(idx, int) for idx in self.g1]):
                    self.results.pairs = list(map(list, [x for x in itertools.product(self.g1, self.g2) if x[0] != x[1]]))
                    # Remove any duplicates if at_g1 and at_g2 are the same
                    for pair1 in self.results.pairs:
                        for pair2 in self.results.pairs:
                            if pair1[0] == pair2[1] and pair1[1] == pair2[0]:
                                self.results.pairs.remove(pair2)
                    self.results.pairs = np.asarray(self.results.pairs)
                else:
                    raise Exception("at_g1 contains lists, however com is set to False and no angle criterion is provided! Please set 'com = True' if you want to calculate the distance between centre of masses of molecules, please specify angle criterion if you want to filter pairs using angle criterion. Otherwise, please use a list of int as ag_g1.")
            # The case when ag_g1 is a H-bond donor list
            else:
                # The case of one donor in at_g1
                if isinstance(self.g1[0], int):
                    try:
                        self.results.pairs = list(set(itertools.product([self.g1[0]], self.g1[1], self.g2)))
                    except TypeError:
                        self.results.pairs = list(set(itertools.product([self.g1[0]], self.g1[1], [self.g2])))
                # The case there are many donors in at_g1
                else:
                    self.results.pairs = []
                    for donor in self.g1:
                        for acceptor in self.g2:
                            if donor[0] == acceptor:
                                symbol = self.traj.frames[0][acceptor].symbs
                                print(f"Donor and acceptor is the same atom: {symbol} {acceptor}, skipped.")
                                continue
                            else:
                                [self.results.pairs.append(x) for x in list(map(list, itertools.product([donor[0]], donor[1], [acceptor])))]
                    self.results.pairs = np.asarray(self.results.pairs)
        
    def pair_filter(self, mic = True, filename = None):
        """
        Calculates the distance (and angle) and compares with the criteria to filter the pairs.
        To save the filtered pairs, use save = True.
        """
        self.results.dimers = []
        # com = self.com
        if self.com:
            for i, frame in enumerate(self.traj.frames):
                if mic:
                    dists = [mic_dist(frame[pair[0].tolist()].calc_com(), frame[pair[1].tolist()].calc_com(), cell = self.traj.cell) for pair in self.results.pairs]
                else:
                    dists = np.linalg.norm(np.asarray([frame[pair[1].tolist()].calc_com() - frame[pair[0].tolist()].calc_com() for pair in self.results.pairs]), axis = 1)
                idx = dists <= self.criteria.distance
                self.results.dimers.append(
                    {
                        "frame": i,
                        "pairs": self.results.pairs[idx].squeeze().tolist(),
                        "dists": dists[idx].tolist(),
                    }
                )
                update_progress(i/self.traj.nframes)
                
        else:
            if isinstance(self.criteria.distance, str):
                temp = re.compile("([\d\.\d]+)([a-zA-Z\_*]+)")
                res = temp.match(self.criteria.distance).groups()
                if res[1] in ["vdw", "vdw_sum", "vdw_r"]:
                    self.criteria.distance = np.asarray([float(res[0]) * (data.vdW_R[self.traj.frames[0][pair[0].tolist()].symbs[0]] + data.vdW_R[self.traj.frames[0][pair[1].tolist()].symbs[0]]) for pair in self.results.pairs])
                elif res[1] in ["coval", "coval_sum", "coval_r"]:
                    self.criteria.distance = np.asarray([float(res[0]) * (data.val_R[self.traj.frames[0][pair[0].tolist()].symbs[0]] + data.val_R[self.traj.frames[0][pair[1].tolist()].symbs[0]]) for pair in self.results.pairs])
                else:
                    raise Exception("Unrecognised string specifier for distance criterion! Please use float + 'vdw_sum' or float + 'coval_sum' as distance criterion, eg '1.05vdw_sum' or '1.5coval_sum'.")
                # Dimers based on distance only
            if self.criteria.angle is None:
                for i, frame in enumerate(self.traj.frames):
                    if mic:
                        dists = np.asarray([frame.dist(pair[0], pair[1], mic = True) for pair in self.results.pairs])
                    else:
                        dists = np.linalg.norm(frame[self.results.pairs[:,1].tolist()].pos - frame[self.results.pairs[:,0].tolist()].pos, axis = 1)
                    idx = dists <= self.criteria.distance
                    self.results.dimers.append(
                        {
                            "frame": i,
                            "pairs": self.results.pairs[idx].squeeze().tolist(),
                            "dists": dists[idx].tolist(),
                        }
                    )
                    update_progress(i/self.traj.nframes)
            # Complicated case, using distance and angle criteria
            else:
                dt = np.dtype([('distance', np.float64), ('angle', np.float64)])
                if isinstance(self.criteria.angle, str):
                    temp = re.compile("([\d\.\d]+)([a-zA-Z]+)")
                    res = temp.match(self.criteria.angle).groups()
                    
                    for i, frame in enumerate(self.traj.frames):
                        if mic:
                            dist_ang = np.asarray([(frame.dist(pair[0], pair[2], mic = True), frame.angle(*pair.tolist(), mic = True)) for pair in self.results.pairs], dtype = dt)
                            d_h_dists = np.asarray([frame.dist(*pair[:-1].tolist(), mic = True) for pair in self.results.pairs])
                            # Add a fluctuation to the angle around the specified angle
                            self.criteria.angle =float(res[0]) * d_h_dists / (1.05 * np.array([data.val_R[self.traj.frames[0][pair[0].tolist()].symbs[0]] + data.val_R[self.traj.frames[0][pair[1].tolist()].symbs[0]] for pair in self.results.pairs]))
                        else:
                            dist_ang = np.asarray([(frame.dist(pair[0], pair[2], mic = False), frame.angle(*pair.tolist(), mic = False)) for pair in self.results.pairs], dtype = dt)
                            d_h_dists = np.asarray([frame.dist(*pair[:-1].tolist(), mic = False) for pair in self.results.pairs])
                            # Add a fluctuation to the angle around the specified angle
                            self.criteria.angle =float(res[0]) * d_h_dists / (1.05 * np.array([data.val_R[self.traj.frames[0][pair[0].tolist()].symbs[0]] + data.val_R[self.traj.frames[0][pair[1].tolist()].symbs[0]] for pair in self.results.pairs]))
                        idx = np.where((dist_ang["distance"] <= self.criteria.distance) & (dist_ang["angle"] >= self.criteria.angle))
                        self.results.dimers.append(
                            {
                                "frame": i,
                                "d_h_a_pairs": self.results.pairs[idx].squeeze().tolist(),
                                "d_a_dists": dist_ang["distance"][idx].tolist(),
                                "d_h_a_angles": dist_ang["angle"][idx].tolist(),
                            }
                        )
                        update_progress(i/self.traj.nframes)
                else:
                    for i, frame in enumerate(self.traj.frames):
                        if mic:
                            dist_ang = np.asarray([(frame.dist(pair[0], pair[2], mic = True), frame.angle(*pair.tolist(), mic = True)) for pair in self.results.pairs], dtype = dt)
                            d_h_dists = np.asarray([frame.dist(*pair[:-1].tolist(), mic = True) for pair in self.results.pairs])
                        else:
                            dist_ang = np.asarray([(frame.dist(pair[0], pair[2], mic = False), frame.angle(*pair.tolist(), mic = False)) for pair in self.results.pairs], dtype = dt)
                            d_h_dists = np.asarray([frame.dist(*pair[:-1].tolist(), mic = False) for pair in self.results.pairs])
                        idx = np.where((dist_ang["distance"] <= self.criteria.distance) & (dist_ang["angle"] >= self.criteria.angle))
                        self.results.dimers.append(
                            {
                                "frame": i,
                                "d_h_a_pairs": self.results.pairs[idx].squeeze().tolist(),
                                "d_a_dists": dist_ang["distance"][idx].tolist(),
                                "d_h_a_angles": dist_ang["angle"][idx].tolist(),
                            }
                        )
                        update_progress(i/self.traj.nframes)
                        
        if filename is not None:
            import json
            with open(filename, "w") as json_file:
                json.dump(self.results.dimers, json_file, indent=4, separators=(", ", ": "))
            json_file.close()
        update_progress(1)


class DLDDF(object):
    """Dimer Lifetime-Displacement Distribution Function
    Takes the results.dimers of a dimers object and calculate the dimer lifetime-displacement distribution function
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Arguments
    at_g1 : list of int, selected atoms to be calculated
    at_g2 : list of int, selected atoms to be calculated
    tau : the time interval
    nbins : int (optional), number of bins (resolution) in the histogram
    criteria : list of geometry criteria for the specified atoms to be regarded as pairs
    com : wether to use the centre of mass of the atom groups in the distance calculation. Set com = True if the atoms groups are two lists of molecules.
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Returns a dataclass containing the analysis results
    x :  np.array, radii over which g(r) is computed
    rdf : 
    plot : 
    t : 
    pairs : 
    scaler : np.array distances
    com_plot :
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Note:
    Periodic boundary conditions not implemented in VHCF, due to the fact we cares the absolute displacement in this analysis.
    """
    def __init__(self, traj, dimer_settings = None, nbins = 100, range = np.array([(0.0, 2.0), (0.0, 5.0)])):
        self.traj = traj
        self.dimer_settings = dimer_settings
        settings_dc = make_dataclass("settings", "bins, range")
        self.settings = settings_dc(nbins, range)
        results = make_dataclass("Results", "dimers r_edges tau_edges r tau dists taus dlddf plot")
        self.results = results
        self.results.dlddf, self.results.r_edges, self.results.tau_edges = np.histogram2d([-1], [-1], **asdict(self.settings))
        self.results.r = (self.results.r_edges[:-1] + self.results.r_edges[1:])/2
        self.results.tau = (self.results.tau_edges[:-1] + self.results.tau_edges[1:])/2
        self.results.dists = []
        self.results.taus = []
        
    def dimer_res(self, mic = True, filename = None):
        dimers_data = dimers(self.traj, **self.dimer_settings)
        dimers_data.make_pairs()
        dimers_data.pair_filter(mic = mic, filename = filename)
        self.results.dimers = dimers_data.results.dimers
    
    def calculate(self, nframes = 100, mic = True, plot = False, bins = 100, levels = 25, xlims = [0, 0.2], ylims = [0, 0.8], filename = None, **kwargs):
        """
        Calculates the dlddf. Set skip = 1 whenever possible to ensure best accuracy."""
        # Use the dimers object to resolve dimers
        key = "pairs"
        if key not in list(self.results.dimers[0].keys()):
            key = "d_h_a_pairs"
        
        dimer_chunks = to_sublists(self.results.dimers, nframes)
        
        for i, dimer_chunk in enumerate(dimer_chunks):
            select = dimer_chunk
            pairs_0 = select[0][key]
            # Real copy the pairs_0 list so that the pairs_0 do not change during mutation
            pairs_ref = pairs_0[:]
            frame_idx_0 = select[0]["frame"]

            # Pass if the first frame in selected dimer_chunk is empty
            if pairs_0 == []:
                # print("Empty pool, pass")
                continue

            else:
                # Iterate through the frames to identify the lifetime and displacement of the pair
                for j, dimer_list in enumerate(select[1:]):
                    pairs = dimer_list[key]
                    frame_idx = dimer_list["frame"]
                    # No pairs in current frame, all pairs are dead
                    if pairs == []:
                        dead_pair = pairs_0[:]
                        # print("All pairs dead!")

                        # Only one pair in the pool
                        if all(isinstance(x, int) for x in pairs_0):
                            if pairs_0 == []:
                                break
                            else:
                            # Calculate the displacement and lifetime
                                if mic:
                                    dist = mic_dist(self.traj.frames[frame_idx_0][dead_pair].calc_com(), self.traj.frames[frame_idx][dead_pair].calc_com(), cell = self.traj.cell)
                                else:
                                    dist = np.linalg.norm(self.traj.frames[frame_idx][dead_pair].calc_com() - self.traj.frames[frame_idx_0][dead_pair].calc_com())
                                lifetime = (j + 1) * self.traj.timestep/1000
                                self.results.dists.append(dist)
                                self.results.taus.append(lifetime)

                        # More than one pairs in the pool
                        else:
                            for dead_pair in dead_pairs:
                                if mic:
                                    dist = mic_dist(self.traj.frames[frame_idx_0][dead_pair].calc_com(), self.traj.frames[frame_idx][dead_pair].calc_com(), cell = self.traj.cell)
                                else:
                                    dist = np.linalg.norm(self.traj.frames[frame_idx][dead_pair].calc_com() - self.traj.frames[frame_idx_0][dead_pair].calc_com())
                                lifetime = (j + 1) * self.traj.timestep/1000
                                self.results.dists.append(dist)
                                self.results.taus.append(lifetime) 

                        break
                    # if 
                    # If there is only one pair in current frame
                    elif all(isinstance(x, int) for x in pairs):
                        if pairs_0 == []:
                            break
                        if all(isinstance(x, int) for x in pairs_0):
                            if pairs_0 == pairs:
                                pass
                        else:
                            dead_pairs = [x for x in pairs_0 if x != pairs]

                            # Calculate the displacement and lifetime
                            for dead_pair in dead_pairs:
                                if mic:
                                    dist = mic_dist(self.traj.frames[frame_idx_0][dead_pair].calc_com(), self.traj.frames[frame_idx][dead_pair].calc_com(), cell = self.traj.cell)
                                else:
                                    dist = np.linalg.norm(self.traj.frames[frame_idx][dead_pair].calc_com() - self.traj.frames[frame_idx_0][dead_pair].calc_com())
                                lifetime = (j + 1) * self.traj.timestep/1000
                                self.results.dists.append(dist)
                                self.results.taus.append(lifetime)

                            [pairs_ref.remove(x) for x in dead_pairs]
                            pairs_0 = pairs_ref[:]

                    # If there are more than one pair in current frame, the pairs will be list of lists
                    else:
                        # Only one pair in the pool
                        if all(isinstance(x, int) for x in pairs_0):
                            if pairs_0 == []:
                                break
                            # Dimmer dead
                            elif pairs_0 not in pairs:
                                # print(f"Dimer {pairs_0} dead, removed from pool.")
                                dead_pair = pairs_0[:]
                                # Calculate the dislpacement and lifetime
                                if mic:
                                    dist = mic_dist(self.traj.frames[frame_idx_0][dead_pair].calc_com(), self.traj.frames[frame_idx][dead_pair].calc_com(), cell = self.traj.cell)
                                else:
                                    dist = np.linalg.norm(self.traj.frames[frame_idx][dead_pair].calc_com() - self.traj.frames[frame_idx_0][dead_pair].calc_com())
                                lifetime = (j + 1) * self.traj.timestep/1000
                                self.results.dists.append(dist)
                                self.results.taus.append(lifetime)

                                pairs_ref = []
                                pairs_0 = pairs_ref[:]

                                break
                            # Dimmer survives    
                            else:
                                pass

                        # More than one pairs exist in the pool
                        else:
                            for item in pairs_0:
                                if item not in pairs:

                                    dead_pair = list(deepflatten(item))
                                    # Calculate the dislpacement and lifetime
                                    if mic:
                                        dist = mic_dist(self.traj.frames[frame_idx_0][dead_pair].calc_com(), self.traj.frames[frame_idx][dead_pair].calc_com(), cell = self.traj.cell)
                                    else:
                                        dist = np.linalg.norm(self.traj.frames[frame_idx][dead_pair].calc_com() - self.traj.frames[frame_idx_0][dead_pair].calc_com())
                                    lifetime = (j + 1) * self.traj.timestep/1000
                                    self.results.dists.append(dist)
                                    self.results.taus.append(lifetime)
                                    try:
                                        pairs_ref.remove(item)
                                        pairs_0 = pairs_ref[:]

                                    except ValueError:
                                        print(item)
                                        print(pairs_0)
                                        print("Value error!")
                                        pass
            update_progress(i/len(dimer_chunks))
        
        x = np.asarray(self.results.taus)
        y = np.asarray(self.results.dists)

        H, xedges, yedges = np.histogram2d(x, y, **asdict(self.settings))

        self.results.dlddf = H = H.T # Transpose so that the plot is correct
        
        if plot:

            fig = plt.figure(figsize=(4.8, 4.0))
            ax  = fig.add_axes([0.16, 0.16, 0.82, 0.75])
            
            xcenters = (xedges[:-1] + xedges[1:])/2
            ycenters = (yedges[:-1] + yedges[1:])/2

            X, Y = np.meshgrid(xcenters, ycenters)

            CS = ax.contourf(X, Y, H, cmap = color_ramp, vmin = 0.01 * H.max(), levels = levels, **kwargs)
            ax.contour(X, Y, H, colors ="k", levels = 25, linewidths = 0.1)

            ax.set_xlabel(r"Lifetime (ps)")
            ax.set_ylabel(r"Displacement (${\AA}$)")
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)

            plt.colorbar(CS, pad=0.035, location = "right", ax = ax)
            
            if filename is not None:
                plt.savefig(filename, dpi = 600)

            plt.show()
            
        update_progress(1)
                
        return self.results
