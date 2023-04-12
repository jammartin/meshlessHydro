#!/usr/bin/env python3

import argparse
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

DIM = 3

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Create an initial condition HDF5 file for a 3D fluid-block  test case.")
    parser.add_argument("--numParticles", "-N", metavar="int", type=int, help="number of particles", required=True)
    parser.add_argument("--twoDimensional", "-t", action="store_true", help="sample 2D 'block'")
    
    args = parser.parse_args()

    N = args.numParticles
    
    outH5 = h5.File("fbN{}.h5".format(N), "w")
    print("Generating fluid block initial conditions with", N, "particles ...")

    if args.twoDimensional:
        DIM = 2
    
    pos = np.empty((N, DIM))
    xv = np.linspace(-.5, .5, round(N**(1./DIM)), endpoint=False)
    yv = np.linspace(-.5, .5, round(N**(1./DIM)), endpoint=False)

    if not args.twoDimensional:
        zv = np.linspace(-.5, .5, round(N**(1./DIM)), endpoint=False)
    
    i = 0
    for x in xv:
        for y in yv:
            if not args.twoDimensional:
                for z in zv:
                    pos[i,0] = x
                    pos[i,1] = y
                    pos[i,2] = z
                    #print("Particle i =", i, "[", x, ",", y, , ",", z, "]")
                    i += 1
            else:
                pos[i, 0] = x
                pos[i, 1] = y
                i += 1
                                    
    # set velocities
    vel = np.zeros(pos.shape)
    # set densities
    rho = np.ones(pos.shape[0])
    # create material ID
    matId = np.zeros(len(rho), dtype=np.int8)
    # volume is 1
    m = rho/N
    # create specific internal energy
    u = np.ones(rho.shape)
    
    outH5.create_dataset("x", data=pos) 
    outH5.create_dataset("v", data=vel)
    outH5.create_dataset("m", data=m)
    outH5.create_dataset("u", data=u)
    outH5.create_dataset("materialId", data=matId)
    
    print("... done. Output written to", outH5.filename)
    outH5.close()
    
