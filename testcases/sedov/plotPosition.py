#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
import os
import re
import argparse
import glob


" based on https://github.com/jammartin/ParaLoBstar/blob/main/tools/conservation/main.py"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot position of particles.")
    parser.add_argument("--data", "-d", metavar="str", type=str, help="input directory",
                        nargs="?", default="output")
    parser.add_argument("--output", "-o", metavar="str", type=str, help="output directory",
                        nargs="?", default="output")

    args = parser.parse_args()

    time = []
    positions = []
    densities = []

for h5file in sorted(glob.glob(os.path.join(args.data, "*.h5")), key=os.path.basename):
    print("Processing ", h5file, " ...")
    data = h5py.File(h5file, 'r')
    time.append(re.findall(r'(\d+)*.h5', h5file))

    print("...reading positions...")
    positions.append(np.array(data["x"][:]))
    densities.append(np.array(data["rho"][:]))
    

print("...done.")

    
numPlots = len(positions)

for i in range(numPlots):
    print("...plotting timestep figure {} / {}....".format(i, numPlots))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    r = positions[i]
    rho = densities[i]
    img = ax.scatter(r[:, 0], r[:, 1], r[:, 2], c=rho)
    fig.colorbar(img, ax=ax)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    #to compare different timesteps
    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)
    ax.set_title('Timestep: {}'.format(time[i][0]))
    
    plt.savefig("{}/Position_at_Timestep{}.png".format(args.output, time[i][0]))
    #plt.show()
