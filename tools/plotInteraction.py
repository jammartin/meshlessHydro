#!/usr/bin/env python3

import argparse
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot interaction geometry for the Riemann problem i->j.")
    parser.add_argument("--fileBasePath", "-f", metavar="string", type=str, help="input file stem", required=True)
    parser.add_argument("--iParticle", "-i", metavar="int", type=int, help="interaction between particle i ...", default=0)
    parser.add_argument("--jParticle", "-j", metavar="int", type=int, help="... and particle j", default=1)

    args = parser.parse_args()
    i = args.iParticle
    j = args.jParticle
    
    plt.rc('text', usetex=True)

    data = h5.File(args.fileBasePath + ".h5", 'r')
    dataNNL = h5.File(args.fileBasePath + "NNL.h5", 'r')

    # retrieve data of particle i from file
    xi = data['x'][i]
    vi = data['v'][i]

    # retrieve data of particle j from file
    xj = data['x'][j]
    vj = data['v'][j]
    
    # finding index ij in NNL data of neighbor j
    nnl = dataNNL["nnl" + str(i)][:]
    ij = np.where(nnl == j)[0][0]

    Aij = dataNNL["Aij" + str(i)][ij]
    vFrame = dataNNL["vFrame" + str(i)][ij]

    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    ax.scatter(xi[0], xi[1], s=50., marker='x', color='b')
    ax.scatter(xj[0], xj[1], s=50., marker='x', color='r')
    ax.quiver(xi[0], xi[1], vi[0], vi[1], angles='xy', scale_units='xy', scale=100.)
    ax.quiver(xj[0], xj[1], vj[0], vj[1], angles='xy', scale_units='xy', scale=100.)

    
    # first order quadrature point
    xij = (xi + xj)/2.
    ax.scatter(xij[0], xij[1], s=50., marker='x', color='k')
    ax.quiver(xij[0], xij[1], vFrame[0], vFrame[1], angles='xy', scale_units='xy', scale=100., color='r')

    ax.quiver(xij[0], xij[1], Aij[0], Aij[1], angles='xy', scale_units='xy', scale=1., color='b')
    
    plt.show()
    
    
    
    
    
    
