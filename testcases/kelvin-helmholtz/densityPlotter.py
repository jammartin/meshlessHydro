#!/usr/bin/env python3

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#MAX_NUM_INTERACTIONS = 1000

def createPlot(h5File, outDir, plotGrad, plotVel, iNNL):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]

    time = data["time"][()][0]
    
    rho = data["rho"][()]
    #P = data["P"][()]
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(7,6), dpi=200)
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=500.) # good for ~100 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=200.) # good for ~400 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=100.) # good for ~900 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=20.) # good for 10**4 particles
    rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=5.) # good for 128**2 particles
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=200.) # good for ~400 particles
    
    if "Ghosts" not in str(h5File):
        #ax.set_xlim((-.75, .75))
        #ax.set_ylim((-.75, .75))
        #ax.set_xlim((-.5, .5))
        #ax.set_ylim((-.5, .5))
        ax.set_xlim((0., 1.))
        ax.set_ylim((0., 1.))

        
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


    plt.title(r"Color coded density $\varrho$ at $t = " + str(time) + r"$")
    #plt.title(r"Color coded pressure $P$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(rhoPlt, cax=cax) #, orientation='horizontal')
    #fig.colorbar(PPlt, ax=ax)
    
    ax.set_aspect('equal')
    
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
    #uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=100.) # good for ~900 particles
    #uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=200.) # good for ~400 particles
    uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=10.) # good for 10**4 particles
    
    if "Ghosts" not in str(h5File):
        ax.set_xlim((-.5, .5))
        ax.set_ylim((-.5, .5))
        #ax.set_xlim((-.6, .6))
        #ax.set_ylim((-.6, .6))

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
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=200.) # good for ~400 particles
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=100.) # good for ~900 particles
    PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=10.) # good for 10**4 particles
    
    if "Ghosts" not in str(h5File):
        ax.set_xlim((-.5, .5))
        ax.set_ylim((-.5, .5))
        #ax.set_xlim((-.6, .6))
        #ax.set_ylim((-.6, .6))
       
    fig.colorbar(PPlt, ax=ax)
    plt.title(r"Color coded pressure $P$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/P" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/P" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def createNoiPlot(h5File, outDir):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    
    noi = data["noi"][()] 
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    noiPlt = ax.scatter(pos[:,0], pos[:,1], c=noi, s=5.) # good for 128**2 particles
    
    if "Ghosts" not in str(h5File):
        ax.set_xlim((-.5, .5))
        ax.set_ylim((-.5, .5))
        #ax.set_xlim((-.6, .6))
        #ax.set_ylim((-.6, .6))
       
    fig.colorbar(noiPlt, ax=ax)
    plt.title(r"Color coded number of interactions")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/noi" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/noi" + pathlib.Path(h5File).stem + ".png")
    plt.close()

   
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot density of results from Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str, help="output directory for generated plots", default="output")
    parser.add_argument("--plotGradient", "-g", action="store_true", help="plot density gradients")
    parser.add_argument("--plotGhosts", "-G", action="store_true", help="also plot ghost cells in an extra file")
    parser.add_argument("--pressure", "-P", action="store_true", help="plot pressure instead of density")
    parser.add_argument("--energy", "-u", action="store_true", help="plot internal energy instead of density")
    parser.add_argument("--noi", "-n", action="store_true", help="plot number of interactions instead of density")
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
                elif args.noi:
                    if args.plotGradient or args.iNNL > -1:
                        print("WARNING: command line arguments '--plotGradient' and '--iNNL' are ignored when plotting energy.")
                    createNoiPlot(h5File, args.outDir)
                else:
                    createPlot(h5File, args.outDir, args.plotGradient, args.plotVelocity, args.iNNL)
            

    print("... done.")
    
