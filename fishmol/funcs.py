"""A range of functions to analyse the trajectory object"""

# RDF
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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