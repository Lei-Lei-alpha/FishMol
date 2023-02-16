"""A range of functions to analyse the trajectory object"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from colour import Color
from matplotlib.colors import LinearSegmentedColormap
from fishmol.utils import to_sublists



ramp_colors = ["#ffffff", "#9ecae1", "#2166ac", "#1a9850", "#ffff33", "#b2182b", "#67000d"]
color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )


# RDF
class RDF(object):
    """
    Pair distribution function:
    Calculated by histogramming distances between all particles in `g1` and `g2` while taking
    periodic boundary conditions into account via the minimum image
    convention.

    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Arguments
    at_g1 : list of int, selected atoms to be calculated
    at_g2 : list of int, selected atoms to be calculated
    nbins : int (optional), number of bins (resolution) in the histogram
    range : tuple or list (optional), the size of the RDF
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Returns a dataclass containing the analysis results
    x :  np.array, radii over which g(r) is computed
    rdf : 
    plot : 
    t : 
    pairs : 
    scaler : np.array distances
    com_plot : 
    """
    def __init__(self, traj, at_g1, at_g2, nbins = 200, range = (0.0, 15.0)):
        self.g1 = at_g1
        self.g2 = at_g2
        self.traj = traj
        settings = make_dataclass("Settings", "bins range")
        self.settings = settings(nbins, range)
        results = make_dataclass("Results", "edges count x rdf plot t pairs scaler com_plot")
        self.results = results
        self.results.count, self.results.edges = np.histogram([-1], **asdict(self.settings))
        self.results.x = (self.results.edges[:-1] + self.results.edges[1:])/2
        self.results.t = np.linspace(0, self.traj.timestep * (self.traj.nframes - 1)/1000, self.traj.nframes)
        # Need to know average volume
        self.volume = np.zeros(self.settings.bins) + np.asarray(self.traj.cell[0]).dot(np.cross(np.asarray(self.traj.cell[1]), np.asarray(self.traj.cell[2])))
        
    def calculate(self, plot = False, com_plot = False, **kwargs):
        """
        Calculate the RDF of two atom groups by calling this function.
         
        """
        dists = np.asarray([self.traj.frames[i].dists(self.g1, self.g2, cutoff = self.settings.range[1], mic = True)[1] for i in range(self.traj.nframes)])
        count = np.histogram(dists, **asdict(self.settings))[0]
        self.results.count = count
        # Use the volume of the simulation box
        
        # Number of each selection
        nA = len(self.g1)
        nB = len(self.g2)
        N = nA * nB
            
        # Volume in each radial shell
        vols = np.power(self.results.edges, 3)
        vol = 4/3 * np.pi * np.diff(vols)
        # Average number density
        density = N / self.volume # number of particles per volume
        # Save pairs as label
        self.results.pairs = list(itertools.product(self.g1, self.g2))
        # Distances for the temporal plot
        self.results.scaler = dists
        # rdf
        self.results.rdf = np.asarray(self.results.count / (density * vol * self.traj.nframes))
        
        if plot:
            fig, ax = plt.subplots(figsize = (4.2, 3.6))
            ax.plot(self.results.x, self.results.rdf, **kwargs)
            ax.set_xlabel(r"$r$ ($\AA$)")
            ax.set_ylabel(r"$g(r)$")
            plt.show()
            self.results.plot = fig
        
        if com_plot:
            fig = plt.figure(figsize=(4.8,3.6))
            ax  = fig.add_axes([0.20, 0.16, 0.685, 0.75])
            for i in range(self.results.scaler.shape[1]):
                ax.plot(self.results.t, self.results.scaler[:,i], label = f"Pair {self.results.pairs[i]}", **kwargs)

            ax.set_xlabel(r"$t$ (ps)")
            ax.set_ylabel(r"$r$ ($\mathrm{\AA}$)")
            plt.legend(frameon = False, ncol = 2, bbox_to_anchor=(0.4, 0.9, 0.4, 0.5), loc='center', fontsize = "small")
            divider = make_axes_locatable(ax)
            ax_histy = divider.append_axes("right", 0.5, pad=0.05, sharey=ax)
            ax_histy.yaxis.tick_right()
            ax_histy.plot(self.results.rdf, self.results.x)
            plt.show()
            self.results.com_plot = fig
            
        return self.results

# ADF
class ADF(object):
    """Angular distribution function
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Arguments
    at_g1 : AtomGroup
      First AtomGroup
    at_g2 : AtomGroup
      Second AtomGroup
    nbins : int (optional)
          Number of bins in the histogram
    range : tuple or list (optional)
          The size of the RDF
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Returns
    g_r : (n_radii) np.array
        radial distribution function values g(r).
    radii : (n_radii) np.array
        radii over which g(r) is computed
    """
    def __init__(self, traj, at_g1, at_g2, at_g3, nbins = 200, range = (0.0, 180.0)):
        self.g1 = at_g1
        self.g2 = at_g2
        self.g3 = at_g3
        self.traj = traj
        settings = make_dataclass("Settings", "bins range")
        self.settings = settings(nbins, range)
        count, edges = np.histogram([-1], **asdict(self.settings))
        count = count.astype(np.float64)
        count *= 0.0
        results = make_dataclass("Results", "edges count x adf plot t pairs scaler com_plot")
        self.results = results
        self.results.edges = edges
        self.results.x = (edges[:-1] + edges[1:])/2
        self.results.count = count
        self.results.t = np.linspace(0, self.traj.timestep * (self.traj.nframes - 1)/1000, self.traj.nframes)
        
    def calculate(self, cone_correction = True, plot = False, com_plot = False, **kwargs):
        angles = []
        for frame in self.traj:
            pairs, angle = frame.angles(self.g1, self.g2, self.g3, mic = True)
            count = np.histogram(angle, **asdict(self.settings))[0]
            self.results.count += count
            angles.append(angle)
        
        self.results.pairs = pairs
        self.results.scaler = np.asarray(angles)
        # Number of each selection
        NA = len(self.g1)
        NB = len(self.g2)
        NC = len(self.g3)
        
        # Calculate ADF
        if cone_correction:
            self.results.adf = np.asarray(self.results.count / (NA * NB * NC * np.sin(self.results.x * np.pi / 180))) # cone correction, convert degree to radians
        else:
            self.results.adf = np.asarray(self.results.count / (NA * NB * NC))
        
        if plot:
            fig, ax = plt.subplots(figsize = (4.2, 3.6))
            ax.plot(self.results.x, self.results.adf, **kwargs)
            ax.set_xlabel(r"$\alpha$ ($^{\circ}$)")
            ax.set_ylabel(r"$g(\alpha)$")
            plt.show()
            self.results.plot = fig
        
        if com_plot:
            fig = plt.figure(figsize=(4.8,3.6))
            ax  = fig.add_axes([0.20, 0.16, 0.685, 0.75])
            for i in range(self.results.scaler.shape[1]):
                ax.plot(self.results.t, self.results.scaler[:,i], label = f"Pair {self.results.pairs[i]}", **kwargs)

            ax.set_xlabel(r"$t$ (ps)")
            ax.set_ylabel(r"$\alpha$ ($^{\circ}$)")
            plt.legend(frameon = False, ncol = 2, bbox_to_anchor=(0.4, 0.9, 0.4, 0.5), loc='center', fontsize = "small")
            divider = make_axes_locatable(ax)
            ax_histy = divider.append_axes("right", 0.5, pad=0.05, sharey=ax)
            ax_histy.yaxis.tick_right()
            ax_histy.plot(self.results.adf, self.results.x)
            plt.show()
            self.results.com_plot = fig
            
        return self.results

# DDF
class DDF(object):
    """Dihedral distribution function
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Arguments
    at_g1 : AtomGroup
      First AtomGroup
    at_g2 : AtomGroup
      Second AtomGroup
    nbins : int (optional)
          Number of bins in the histogram
    range : tuple or list (optional)
          The size of the RDF
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Returns
    g_r : (n_radii) np.array
        radial distribution function values g(r).
    radii : (n_radii) np.array
        radii over which g(r) is computed
    """
    def __init__(self, traj, at_g, nbins = 200, range = (0.0, 180.0)):
        self.g = at_g
        self.traj = traj
        settings = make_dataclass("Settings", "bins range")
        self.settings = settings(nbins, range)
        count, edges = np.histogram([-1], **asdict(self.settings))
        count = count.astype(np.float64)
        count *= 0.0
        results = make_dataclass("Results", "edges count x ddf plot t angles com_plot")
        self.results = results
        self.results.edges = edges
        self.results.x = (edges[:-1] + edges[1:])/2
        self.results.count = count
        self.results.t = np.linspace(0, self.traj.timestep * (self.traj.nframes - 1)/1000, self.traj.nframes)
        
    def calculate(self, plot = False, com_plot = False, **kwargs):
        angles = []
        for frame in self.traj:
            pairs, angle = frame.dihedrals(self.g, mic = True)
            count = np.histogram(angle, **asdict(self.settings))[0]
            self.results.count += count
            angles.append(angle)
        
        self.results.angles = np.asarray(angles)
        
        N = len(pairs)
        
        # Calculate DDF
        self.results.ddf = np.asarray(self.results.count) / N
        
        if plot:
            fig, ax = plt.subplots(figsize = (4.2, 3.6))
            ax.plot(self.results.x, self.results.ddf, **kwargs)
            ax.set_xlabel(r"$\delta$ ($^{\circ}$)")
            ax.set_ylabel(r"$g(\delta)$")
            plt.show()
            self.results.plot = fig
        
        if com_plot:
            fig = plt.figure(figsize=(4.8,3.6))
            ax  = fig.add_axes([0.20, 0.16, 0.685, 0.75])
            for i in range(self.results.angles.shape[1]):
                ax.plot(self.results.t, self.results.angles[:,i], **kwargs)

            ax.set_xlabel(r"$t$ (ps)")
            ax.set_ylabel(r"$\delta$ ($^{\circ}$)")
            divider = make_axes_locatable(ax)
            ax_histy = divider.append_axes("right", 0.5, pad=0.05, sharey=ax)
            ax_histy.yaxis.tick_right()
            ax_histy.plot(self.results.ddf, self.results.x)
            plt.show()
            self.results.com_plot = fig
            
        return self.results


# CDF
class CDF(object):
    """Combined distribution function
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Arguments
    at_g1 : AtomGroup
      First AtomGroup
    at_g2 : AtomGroup
      Second AtomGroup
    nbins : int (optional)
          Number of bins in the histogram
    range : tuple or list (optional)
          The size of the RDF
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Returns
    g_r : (n_radii) np.array
        radial distribution function values g(r).
    radii : (n_radii) np.array
        radii over which g(r) is computed
    """
    def __init__(self, scaler1, scaler2, names = None, range = None):
        self.s1 = scaler1
        self.s2 = scaler2
        self.names = names
        if range is None:
            range = ([0, getattr(self.s1, names[0]).max()],[0, getattr(self.s2, names[1]).max()])
        settings = make_dataclass("Settings", "range")
        self.settings = settings(range)
        results = make_dataclass("Results", "xedges yedges xcenters ycenters x y cdf plot")
        self.results = results
        
    def calculate(self, plot = True, **kwargs):
        self.results.xedges = self.s1.edges
        self.results.yedges = self.s2.edges

        x = self.s1.scaler.flatten()
        y = self.s2.scaler.flatten()

        self.results.cdf, self.results.xedges, self.results.yedges = np.histogram2d(x, y, bins = (self.results.xedges, self.results.yedges), **asdict(self.settings))
        self.results.cdf = self.results.cdf.T # Transpose so that the plot is correct
        
        if plot:
            fig = plt.figure(figsize=(4.8, 4.0))
            ax  = fig.add_axes([0.16, 0.16, 0.82, 0.75])

            X, Y = np.meshgrid(self.s1.x, self.s2.x)
            levels = np.linspace(self.results.cdf.min(), self.results.cdf.max(), 50)
            CS = ax.contourf(X, Y, self.results.cdf, cmap = color_ramp, levels = levels)

            # ax.set_xlabel(r"$r$ ($\AA$)")
            # ax.set_ylabel(r"$\alpha$ ($^{\circ}$)")
            # ax.set_xlim(2,5)

            divider = make_axes_locatable(ax)

            ax_histx = divider.append_axes("top", 0.5, pad=0.05, sharex=ax)
            ax_histx.xaxis.set_tick_params(labelbottom=False)
            ax_histx.yaxis.set_tick_params(labelleft=False)
            ax_histx.xaxis.tick_top()
            ax_histx.plot(self.s1.x, getattr(self.s1, self.names[0]))

            ax_histy = divider.append_axes("right", 0.5, pad=0.05, sharey=ax)
            ax_histy.xaxis.set_tick_params(labelbottom=False)
            ax_histy.yaxis.set_tick_params(labelleft=False)
            ax_histy.yaxis.tick_right()
            ax_histy.plot(getattr(self.s2, self.names[1]), self.s2.x)

            plt.colorbar(CS, pad=0.035, location = "right", ax = ax)
            # plt.savefig("cdf_water14_15_water143.jpg", dpi = 600)
            self.results.plot = fig
            plt.show()
            
        return self.results


# VRD
class VRD(object):
    """Vector Reorientation Dynamics
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Arguments
    spec : np.array or list of two atoms (two atom groups) that determine the vector
    timestep : the timestep of the vector when using np.array input
    traj : the input Trajectory object
    num : int, length of frames to calculate the reorientation dynamics
    skip : int, frames to skip = skip - 1. Determines the sampling step of the Trajectory, eg. skip = 1, use all frames in the Trajectory: skip = 10, skip 9 frames.
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Returns
    sdf : np.array, scaler distribution function values g.
    radii : np.array, radii over which g is computed
    """
    def __init__(self, traj = None, spec = None, timestep = None, num = 2000, sampling = 5, skip = 10):
        results = make_dataclass("Results", "t c_t_mean c_t_error plot")
        self.results = results
        self.traj = traj
        self.spec = spec
        self.num = num
        self.sampling = sampling
        self.skip = skip
        if traj is not None:
            self.traj = traj
            self.t_step = traj.timestep
        elif traj is None:
            if isinstance(self.spec, np.ndarray):
                if timestep is not None:
                    self.t_step = timestep
                else:
                    print("VRD created with array of vectors without specifying timestep, using default timestep of 5 fs, please specify the timestep you are using if this is incorrect.")
                    self.t_step = 5
            else:
                raise Exception("Please either specify\n:(i)a Trajectory object and a list of the indices of two atoms\nor\n(ii) a np.array object of vector data and the timestep of the vector data.")
        else:
            raise Exception("Please either specify\n:(i)a Trajectory object and a list of the indices of two atoms\nor\n(ii) a np.array object of vector data and the timestep of the vector data.")
        self.results.t = np.linspace(0, self.t_step * (self.num - 1) / 1000, num = self.num//self.sampling)
        
    def calculate(self, plot = True, log_scale = False, **kwargs):
        if self.traj is not None:
            frame_chunks = to_sublists(self.traj.frames, self.num)[::self.skip]
            dot_products = []
            for i, frame_chunk in enumerate(frame_chunks):
                select = frame_chunk[::self.sampling]
                # (a * (b)).sum(axis=1)
                dot_products.append([frame.vecs(*self.spec, absolute = False, normalise = True, mic = True).flatten().dot(
                    select[0].vecs(*self.spec, absolute = False, normalise = True, mic = True).flatten().T
                ) for frame in select])
                update_progress(i / len(frame_chunks))
            dot_products = np.asarray(dot_products)
            
            update_progress(1)
            self.results.c_t_mean = dot_products.mean(axis = 0)
            self.results.c_t_error = dot_products.std(axis = 0) / (len(frame_chunks))**0.5
            
        # If the vec is an array of vectors without traj
        else:
            vec_chunks = to_sublists(self.spec, self.num)[::self.skip]
            dot_products = []
            for i, vec_chunk in len(vec_chunks):
                select = vec_chunk[::self.sampling]
                dot_products.append([vecs.dot(vec_chunk[0]) / (np.linalg.norm(vecs) * vecs(vec_chunk[0])) for vecs in vec_chunk])
                update_progress(i / len(frame_chunks))
        
        update_progress(1)
        dot_products = np.asarray(dot_products)
        self.results.c_t_mean = dot_products.mean(axis = 0)
        self.results.c_t_error = dot_products.std(axis = 0) / (len(frame_chunks))**0.5
        
        # delete temp variable to release some memory
        del frame_chunks, dot_products, select
        
        # Plot the results
        
        if plot:
            fig, ax = plt.subplots(figsize = (4.2, 3.6))
            ax.plot(self.results.t, self.results.c_t_mean, **kwargs)
            ax.set_xlabel(r"$t (ps)$")
            ax.set_ylabel(r"$C_t$")
            if log_scale:
                plt.semilogy()
            plt.show()
            self.results.plot = fig
            
        return self.results


