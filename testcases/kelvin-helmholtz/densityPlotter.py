#!/usr/bin/env python3

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

#MAX_NUM_INTERACTIONS = 1000

def createPlot(h5File, outDir, plotGrad, iNNL):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    rho = data["rho"][()]
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=500.) # good for ~100 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=10.) # good for 10**4 particles
    
    # plot gradient
    if plotGrad:
        plotGradient(data["rhoGrad"][:], pos, ax);

    # plot NNL for particle i
    if iNNL > -1 and "Ghosts" not in str(h5File):
        plotNNL(h5File, iNNL, pos, ax)
        
    fig.colorbar(rhoPlt, ax=ax)
    plt.title(r"Color coded density $\rho$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/" + pathlib.Path(h5File).stem + ".png")
    plt.close()
    #plt.show()

def plotGradient(grad, pos, ax):
    ax.quiver(pos[:,0], pos[:,1], grad[:,0], grad[:,1], angles='xy', scale_units='xy', scale=1.)
    #ax.quiver(pos[:,0], pos[:,1], grad[:,0], grad[:,1], angles='xy', scale=.01)
    
    #for i, rhoGrad in enumerate(grad):
    #    if np.linalg.norm(rhoGrad) > .05:
    #        print("rhoGrad @", i, "=", rhoGrad)

def plotNNL(h5File, iNNL, pos, ax):
    data = h5.File(str(h5File).replace(".h5", "NNL.h5"), 'r')
    posNNL = data["nnlPrtcls"+str(iNNL)][:]
    ax.scatter(pos[iNNL,0], pos[iNNL,1], s=50., marker='x', color='b')
    ax.scatter(posNNL[:,0], posNNL[:,1], s=50, marker='x', color='r')

   
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot density of results from Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str, help="output directory for generated plots", default="output")
    parser.add_argument("--plotGradient", "-g", action="store_true")
    parser.add_argument("--iNNL", "-i", metavar="int", type=int, help="plot NNL for particles i", default=-1)

    args = parser.parse_args()

    plt.rc('text', usetex=True)
    
    print("Examining files in", args.simOutputDir, "...")
    
    for h5File in pathlib.Path(args.simOutputDir).glob('*.h5'):
        if "NNL" not in str(h5File):
            print("\t", h5File)
            createPlot(h5File, args.outDir, args.plotGradient, args.iNNL)

    print("... done.")
    
