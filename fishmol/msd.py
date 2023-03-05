import numpy as np
from fishmol.utils import update_progress, vector, Arrow3D

def autocorrFFT(x):
    N=len(x)
    F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   #now we have the autocorrelation in convention B
    n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
    return res/n #this is the autocorrelation in convention A

def msd_fft(r):
    N=len(r)
    D=np.square(r).sum(axis=1) 
    D=np.append(D,0) 
    S2=sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        Q=Q-D[m-1]-D[N-m]
        S1[m]=Q/(N-m)
    return S1-2*S2


def traj_proj(pos, vec):
    pos = np.array(pos)
    vec = np.array(vec)
    return np.dot(pos, vec) / np.linalg.norm(vec)

def msd_1d(r):
    shifts = np.arange(len(r))
    msds = np.zeros(r.shape[0])    

    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs)
        msds[i] = sqdist.mean()

    return msds


def proj_d(df, num, start_idx = 0, n_theta = 120, n_phi = 120, timestep = None, filename = None, plot = False, figname =None):
    """
    Calculates the MSD and D alogn the direction specified by vec.
    """
    # A spherical mesh of vectors
    phi = np.linspace(-np.pi, np.pi, num = n_phi)
    theta = np.linspace(-np.pi/2, np.pi, num = n_theta)

    phi, theta = np.meshgrid(phi, theta)
    
    x = np.sin(phi) * np.sin(theta)
    y = np.sin(phi) * np.cos(theta)
    z = np.cos(phi)
    
    vecs = np.stack((x.flatten(), y.flatten(), z.flatten()), axis = 1) 
    ave_D = np.zeros(len(vecs))
    
    if timestep is None:
        timestep = 5 # defauting timestep of trajectory to 5 fs
    t0 = start_idx * timestep
    t_end = t0 + timestep * (len(df) - 1) 
    t = np.linspace(t0, t_end, num = len(df))
    t = t[1000:] - 5000
    
    for i, vec in enumerate(vecs):
        msd_df = np.zeros((len(df), num))
        for j in range(num):
            temp = np.array(df.iloc[:, 3*j:3*j + 3])
            temp = traj_proj(temp, vec)
            msds = msd_fft(temp.reshape(len(temp), 1))
            msd_df[:, j] = msds

        ave_msd = msd_df.mean(axis=1)[1000:]
        
        linear_model=np.polyfit(t/1000, ave_msd, 1)
        ave_D[i] = linear_model[0]/20000
        update_progress(i / len(vecs))
    update_progress(1)
    return vecs, ave_D