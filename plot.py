import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 13,
#    "font.family": "sans-serif",
#    "font.sans-serif": "Arial",
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "figure.figsize": (4.2,3.6),
    "figure.subplot.left": 0.21,
    "figure.subplot.right": 0.96,
    "figure.subplot.bottom": 0.18,
    "figure.subplot.top": 0.93,
    "legend.frameon": False,
})

color = "#b2182b"

def ave_msd_d(df, filename = None, show = False):
    _, ax = plt.subplots(figsize=(4.6,3.6), constrained_layout=True)

    ax.plot(df["t"]/1000, df["Mean_MSD"], color = "#525252")

    ax.plot(df["t"]/1000, df["Mean_MSD"] - df["MSD_error"], linewidth = 1, color = "#bdbdbd")
    ax.plot(df["t"]/1000, df["Mean_MSD"] + df["MSD_error"], linewidth = 1, color = "#bdbdbd")

    ax.fill_between(df["t"]/1000, df["Mean_MSD"] + df["MSD_error"], df["Mean_MSD"] - df["MSD_error"], color= "#d9d9d9", alpha = 0.2)

    ax1=ax.twinx()

    ax1.plot(df["t"]/1000, df["Mean_D"], color = color)
    ax1.plot(df["t"]/1000, df["Mean_D"] - df["D_error"], linewidth = 1, color = "#d6604d")
    ax1.plot(df["t"]/1000, df["Mean_D"] + df["D_error"], linewidth = 1, color = "#d6604d")

    ax1.fill_between(df["t"]/1000, df["Mean_D"] + df["D_error"], df["Mean_D"] - df["D_error"], color = "#fddbc7", alpha = 0.2)

    ax.set_ylabel('MSD (Å$^2$)')
    ax.set_xlabel('$t$ (ps)')
    ax1.set_ylabel('$D$ (cm$^2$ s$^{-1}$)', color = color)
    ax1.tick_params(axis='y', color = color, labelcolor = color)
    ax1.spines['right'].set_color(color)
    ax1.ticklabel_format(axis='y', style='sci', scilimits=[-4,4], useMathText=True)

    if filename != None:
        plt.savefig(filename,dpi=600)
    if show:
        plt.show()


def msd_d_comp(df, atom_dicts, filename=None, show = False):
    _, ax = plt.subplots(nrows=len(atom_dicts)//2, ncols=2, figsize = (9.0, 3.5 * len(atom_dicts)//2), constrained_layout=True)

    color = "#b2182b"

    for i in range(len(atom_dicts)//2):
        for j in range(2):
            k = 2*i+j+1
            
            ax[i, j].plot(df["t"]/1000, df[f"water{k}_MSD"], color = "#525252")
                    
            globals()[f"ax1_{i}_{j}"] = ax[i, j].twinx()
            globals()[f"ax1_{i}_{j}"].plot(df["t"]/1000, df[f"water{k}_D"], color = color)
            
            x_min, x_max = ax[i, j].get_xlim()
            _, y_max = ax[i, j].get_ylim()
            
            ax[i, j].text((x_min + x_max)/2, 1.05*y_max, f"water{k}", va = "center", ha = "center")
            
            ax[i, j].set_ylabel('MSD (Å$^2$)')
            ax[i, j].set_xlabel('$t$ (ps)')
            globals()[f"ax1_{i}_{j}"].set_ylabel('$D$ (cm$^2$ s$^{-1}$)', color = color)
            globals()[f"ax1_{i}_{j}"].tick_params(axis='y', color = color, labelcolor = color)
            globals()[f"ax1_{i}_{j}"].spines['right'].set_color(color)
            globals()[f"ax1_{i}_{j}"].ticklabel_format(axis='y', style='sci', scilimits=[-4,4], useMathText=True)
    if filename != None:
        plt.savefig(filename, dpi =600)
    if show:
        plt.show()


def plot_comp(df, index, direction, filename = None, dpi=600, show = False):
    _, ax = plt.subplots(nrows=2, ncols=2, figsize = (9.0, 7.0), constrained_layout=True)

    for i in range(2):
        for j in range(2):
            k = 2*i+j
            ax[i, j].plot(df.iloc[999:,0]/1000, df.iloc[999:, 8*(index-1) + 2*k + 1], color = "#525252")

            globals()[f"ax1_{i}_{j}"] = ax[i, j].twinx()
            globals()[f"ax1_{i}_{j}"].plot(df.iloc[999:,0]/1000, df.iloc[999:, 8*(index-1) + 2*k + 2], color = color)
            
            x_min, x_max = ax[i, j].get_xlim()
            
            _, y_max = ax[i, j].get_ylim()
        
            ax[i, j].text((x_min+x_max)/2, 1.05*y_max, direction[k], va = "center", ha = "center")

            ax[i, j].set_ylabel('MSD (Å$^2$)')
            ax[i, j].set_xlabel('$t$ (ps)')
            
            globals()[f"ax1_{i}_{j}"].set_ylabel('$D$ (cm$^2$ s$^{-1}$)', color = color)
            globals()[f"ax1_{i}_{j}"].tick_params(axis='y', color = color, labelcolor = color)
            globals()[f"ax1_{i}_{j}"].spines['right'].set_color(color)
            globals()[f"ax1_{i}_{j}"].ticklabel_format(axis='y', style='sci', scilimits=[-4,4], useMathText=True)
    if filename !=None:
        plt.savefig(filename, dpi=dpi)
    if show:
        plt.show()


def quick_plot(x, y, reorientations=[], flips = [], filename = None, show = False):
    fig, ax = plt.subplots()
    
    ax.plot(x/1000, y, color ="#969696", zorder = 0)
    
    ax.set_xlabel("$t$ (ps)")
    ax.set_ylabel("Rotate angle (°)")
    
    if reorientations != []:
        ax.scatter(np.array(reorientations)[:,0], np.array(reorientations)[:,1], color = "#b2182b", alpha =0.6, linewidths=0, label = "Reorientations", zorder = 1)
        plt.legend()
        
    if flips != []:
        ax.scatter(np.array(flips)[:,0], np.array(flips)[:,1], marker = "s", color = "#2166ac", alpha =0.6, linewidths=0, label = "Flips", zorder = 1)
        plt.legend()
    
    if filename != None:
        plt.savefig(filename, dpi = 600)
    if show:
        plt.show()