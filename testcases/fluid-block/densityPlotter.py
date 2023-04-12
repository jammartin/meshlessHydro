#!/usr/bin/env python3

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

#MAX_NUM_INTERACTIONS = 1000

# plotting parameters, to be adjusted for different number of particles and configurations
LOWER_LIM_X = -1.
UPPER_LIM_X = 1.
LOWER_LIM_Y = -1.
UPPER_LIM_Y = 1.
MARKER_SIZE = 200.


def createPlot(h5File, outDir, plotGrad, plotVel, iNNL):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    
    rho = data["rho"][()]
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=MARKER_SIZE) 
        
    if "Ghosts" not in str(h5File):
        ax.set_xlim((LOWER_LIM_X, UPPER_LIM_X))
        ax.set_ylim((LOWER_LIM_Y, UPPER_LIM_Y))
    
    # Plot gradient
    if plotGrad and not plotVel:
        plotGradient(data["rhoGrad"][:], pos, ax)
    elif not plotGrad and plotVel:
        plotVelocity(data["v"][:], pos, ax)
    elif plotGrad and plotVel:
        print("WARNING: command line arguments '--plotVelocity' and '--plotGradient' are incompatible. - Plotting neither.")

    # plot NNL for particle i
    if iNNL > -1 and "Ghosts" not in str(h5File):
        plotNNL(h5File, iNNL, pos, ax)
        
    fig.colorbar(rhoPlt, ax=ax)
    #fig.colorbar(PPlt, ax=ax)
    plt.title(r"Color coded density $\rho$")
    #plt.title(r"Color coded pressure $P$")
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

def plotVelocity(vel, pos, ax):
    ax.quiver(pos[:,0], pos[:,1], vel[:,0], vel[:,1], angles='xy', scale_units='xy', scale=10.)


def plotNNL(h5File, iNNL, pos, ax):
    data = h5.File(str(h5File).replace(".h5", "NNL.h5"), 'r')
    posNNL = data["nnlPrtcls"+str(iNNL)][:]
    ax.scatter(pos[iNNL,0], pos[iNNL,1], s=50., marker='x', color='b')
    ax.scatter(posNNL[:,0], posNNL[:,1], s=50, marker='x', color='r')

def createEnergyPlot(h5File, outDir):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    
    u = data["u"][()] 
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=MARKER_SIZE) 
    
    if "Ghosts" not in str(h5File):
        ax.set_xlim((LOWER_LIM_X, UPPER_LIM_X))
        ax.set_ylim((LOWER_LIM_Y, UPPER_LIM_Y))
               
    fig.colorbar(uPlt, ax=ax)
    plt.title(r"Color coded internal energy $u$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/u" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/u" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def createPressurePlot(h5File, outDir):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    
    P = data["P"][()] 
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=MARKER_SIZE) 
    
    if "Ghosts" not in str(h5File):
        ax.set_xlim((LOWER_LIM_X, UPPER_LIM_X))
        ax.set_ylim((LOWER_LIM_Y, UPPER_LIM_Y))
        
    fig.colorbar(PPlt, ax=ax)
    plt.title(r"Color coded pressure $P$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/P" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/P" + pathlib.Path(h5File).stem + ".png")
    plt.close()
    
   
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot density of results from Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str, help="output directory for generated plots", default="output")
    parser.add_argument("--plotGradient", "-g", action="store_true", help="plot density gradients")
    parser.add_argument("--plotGhosts", "-G", action="store_true", help="also plot ghost cells in an extra file")
    parser.add_argument("--pressure", "-P", action="store_true", help="plot pressure instead of density")
    parser.add_argument("--energy", "-u", action="store_true", help="plot internal energy instead of density")
    parser.add_argument("--plotVelocity", "-v", action="store_true", help="plot velocity")
    parser.add_argument("--iNNL", "-i", metavar="int", type=int, help="plot NNL for particles i", default=-1)

    args = parser.parse_args()

    plt.rc('text', usetex=True)
    
    print("Examining files in", args.simOutputDir, "...")
    
    for h5File in pathlib.Path(args.simOutputDir).glob('*.h5'):
        if "NNL" not in str(h5File):
            if args.plotGhosts or not "Ghost" in str(h5File): 
                print("\t", h5File)
                if args.pressure:
                    if args.plotGradient or args.iNNL > -1:
                        print("WARNING: command line arguments '--plotGradient' and '--iNNL' are ignored when plotting pressure.")
                    createPressurePlot(h5File, args.outDir)
                elif args.energy:
                    if args.plotGradient or args.iNNL > -1:
                        print("WARNING: command line arguments '--plotGradient' and '--iNNL' are ignored when plotting energy.")
                    createEnergyPlot(h5File, args.outDir)
                else:
                    createPlot(h5File, args.outDir, args.plotGradient, args.plotVelocity, args.iNNL)
            

    print("... done.")
    
