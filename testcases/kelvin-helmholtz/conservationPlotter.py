#!/usr/bin/env python3

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot conserved qunatities over time for the Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--label", "-l", metavar="string", type=str, help="label appended to the output file name", default="Unlabeled")
    parser.add_argument("--MFM", "-M", action="store_true", help="do not plot mass as the error is zero")
    args = parser.parse_args()
    
    plt.rc('text', usetex=True)
    #plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    #\renewcommand\vec[1]{\bm{#1}}')
    
    print("Examining files in", args.simOutputDir, "...")

    mass = []
    energy = []
    momX = []
    momY = []
    time = []

    setReference =  True
    
    for h5File in pathlib.Path(args.simOutputDir).glob('*.h5'):
        data = h5.File(h5File, 'r')
        if setReference:
            refMass = data["totalMass"][0]
            refEnergy = data["energy"][0]
            refMomX = data["xMomentum"][0]
            refMomY = data["yMomentum"][0]
            setReference = False
        time.append(data["time"][0])
        mass.append(abs(data["totalMass"][0]-refMass))
        energy.append(abs(data["energy"][0]-refEnergy))
        momX.append(abs(data["xMomentum"][0]-refMomX))
        momY.append(abs(data["yMomentum"][0]-refMomY))

    print("... plotting ... ")    

    plt.rcParams.update({'font.size': 18})
    
    fig, ax = plt.subplots(figsize=(7,6), dpi=200)

    print("M_tot =",  mass)
    print("pX_tot =", momX)
    print("pY_tot =", momY)
    print("E_tot =", energy)
    
    if not args.MFM:
        plt.plot(time, mass, 'ro', label=r'$\Delta M_\text{tot}$')
    plt.plot(time, momX, 'bv', label=r'$\Delta p_{x, \text{tot}}$')
    plt.plot(time, momY, 'g^', label=r'$\Delta p_{y, \text{tot}}$')
    plt.plot(time, energy, 'kx', label=r'$\Delta E_\text{tot}$')
    plt.yscale('log')
    
    plt.title(r"Absolute numerical error $\Delta_\text{num}$ over time $t$")
    #plt.title(r"Color coded pressure $P$")
    plt.xlabel(r"Time $t$")
    plt.ylabel(r"$\Delta_\text{num}$")

    plt.legend(loc='lower left')
    plt.grid()
    
    plt.tight_layout()
    print("Saving figure to conservation" + args.label + ".png")
    plt.savefig("conservation" + args.label + ".png")
    plt.close()
    
    print("... done.")
    
