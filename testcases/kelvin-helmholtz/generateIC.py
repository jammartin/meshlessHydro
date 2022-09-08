#!/usr/bin/env python3

import argparse
import numpy as np
import h5py as h5

DIM = 2

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Create an initial condition HDF5 file for a 2D Kelvin-Helmholtz test case.")
    parser.add_argument("--numParticles", "-N", metavar="int", type=int, help="number of particles", required=True)
    
    args = parser.parse_args()

    N = args.numParticles

    # initialize default random generator
    rng = np.random.default_rng(6102003)
    
    outH5 = h5.File("kh_N{}.h5".format(N), "w")
    print("Generating Kelvin-Helmholtz initial conditions with", N, "particles ...")

    # randomly generate N points in DIM dimensions
    pos = rng.random(size=(N, DIM))
    

    
