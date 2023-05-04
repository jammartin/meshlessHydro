#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc

#iPsi=0

# cubic spline kernel in 2 dimensions
def cubicSpline(r, h):
    sigma = 10./(7.*np.pi*h**2.) # for dim=2
    q = r/h
    if 0. <= q and q <= 1.:
        return sigma*(1.-3./2.*q**2.*(1.-q/2.))
    elif 1. < q and q < 2.:
        return sigma/4.*(2.-q)**3.
    else:
        return 0.

# x is a 2D vector
def psi_i(i, x, y, h, pos):
    xyv = np.array([x, y]) 
    omega = 0.
    for xy in pos:
        omega = omega + cubicSpline(np.linalg.norm(xyv - xy), h)
        
    return 1./omega * cubicSpline(np.linalg.norm(xyv - pos[i]), h)


if __name__=="__main__":

    # initialize random generator
    rng = np.random.default_rng(6102003)

    # randomly generate 3 points in 2 dimensions
    pos = np.empty((3, 2))
    for i in range(3):
        pos[i] = rng.random(size=2)

    # create periodic boundaries
    posPeriodic = np.empty((9*3, 2))
    posPeriodic[:3] = pos
    posPeriodic[3:6] = pos + np.array([-1., -1.])
    posPeriodic[6:9] = pos + np.array([0., -1.])
    posPeriodic[9:12] = pos + np.array([1., -1.])
    posPeriodic[12:15] = pos + np.array([1., 0.])
    posPeriodic[15:18] = pos + np.array([1., 1.])
    posPeriodic[18:21] = pos + np.array([0., 1.])
    posPeriodic[21:24] = pos + np.array([-1., 1.])
    posPeriodic[24:27] = pos + np.array([-1., 0.])
    
    # create meshgrid for plotting    
    x = np.arange(0., 1., .01)
    y = np.arange(0., 1., .01)
    xv, yv = np.meshgrid(x, y, indexing='ij')

    psiv = np.empty((9*3, 100, 100))

    print("Filling meshgrid ...")
    for iPsi in [0, 1, 2]: #, 8*3, 2*3+2, 4*3+1]:
        for i, xi in enumerate(x):
            for j, yj in enumerate(y):
                psiv[iPsi, i, j] = psi_i(iPsi, xi, yj, .3, posPeriodic)
                for k in range(3,27,3):
                    psiv[iPsi, i, j] += psi_i(k+iPsi, xi, yj, .3, posPeriodic)
    print("... done.")

    print("Plotting ...")
    # create transparent color maps
    c_red = mplc.colorConverter.to_rgba('red')
    c_blue= mplc.colorConverter.to_rgba('blue')
    c_green= mplc.colorConverter.to_rgba('green')
    c_white_trans = mplc.colorConverter.to_rgba('white',alpha = 0.0)

    cmap_r = mplc.LinearSegmentedColormap.from_list('r_cmap',[c_white_trans,c_red],512)
    cmap_b = mplc.LinearSegmentedColormap.from_list('b_cmap',[c_white_trans,c_blue],512)
    cmap_g = mplc.LinearSegmentedColormap.from_list('g_cmap',[c_white_trans,c_green],512)

    plt.rc('text', usetex=True)
    plt.rcParams.update({'font.size': 18})
    
    plt.figure(figsize=(6,6), dpi=200)
    
    # plot real cells
    plt.contourf(xv, yv, psiv[0], levels=100, cmap=cmap_r)
    plt.contourf(xv, yv, psiv[1], levels=100, cmap=cmap_b)
    plt.contourf(xv, yv, psiv[2], levels=100, cmap=cmap_g)

    # plot ghost cells
    #plt.contourf(xv, yv, psiv[8*3], levels=100, cmap=cmap_r)
    #plt.contourf(xv, yv, psiv[2*3+2], levels=100, cmap=cmap_g)
    #plt.contourf(xv, yv, psiv[4*3+1], levels=100, cmap=cmap_b)
    
    # plot particles
    plt.plot(pos[:, 0], pos[:, 1], 'kx')

    plt.xlabel("$x$")
    plt.ylabel("$y$")

    plt.gcf().get_axes()[0].set_aspect('equal')
    
    #plt.show()
    plt.tight_layout()
    plt.savefig("volumePartition.png") 
    print("...done. Saved plot to volumePartition.png")
