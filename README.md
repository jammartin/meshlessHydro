# meshlessHydro
Implementing a meshless scheme for hydrodynamic particle simulations as described implemented in [GIZMO](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html). Below a Kelvin-Helmholtz testcase with 10000 particles simulated with the **meshless finite volume (MFV)** algorithm is shown.

<img src="media/KH_N10000.gif" alt="Kelvin-Helmholtz testcase" width="100%"/>

## Implementation roadmap

- [Volume partition](#volume-partition)

### Gradient estimator
For the method all gradients of a quantity $f$ are calculated consitently with locally-centered least-squares matrix gradient operators. A second-order accurate expression is given by 

$$\left( \nabla f \right)_i^\alpha = 
\sum_j (f_j - f_i)\tilde{\psi}_j^\alpha\left( \vec{x}_i \right)$$

where $\tilde{\psi}_j^\alpha\left( \vec{x}_i \right)$ is somewhat complex helper function.


### Riemann Solver

An exact Riemann solver taken from [this repository](https://github.com/bwvdnbro/python_finite_volume_solver) is used to solve the Riemann problem at the effective faces.


## Resources

- [Paper describing the algorithm](https://arxiv.org/abs/1409.7395)

## Project structure

### File tree
```
├── LICENSE
├── README.md
├── demonstrator
│   ├── Makefile
│   ├── config.info
│   ├── debugLog.txt
│   ├── include
│   │   ├── ConfigParser.h
│   │   ├── Domain.h
│   │   ├── Helper.h
│   │   ├── InitialDistribution.h
│   │   ├── Logger.h
│   │   ├── MeshlessScheme.h
│   │   ├── Particles.h
│   │   ├── Riemann.h
│   │   └── parameter.h
│   ├── log
│   └── src
│       ├── ConfigParser.cpp
│       ├── Domain.cpp
│       ├── Helper.cpp
│       ├── InitialDistribution.cpp
│       ├── Logger.cpp
│       ├── MeshlessScheme.cpp
│       ├── Particles.cpp
│       ├── Riemann.cpp
│       └── main.cpp
├── media
│   ├── KH_N10000.gif
│   └── volumePartition.png
├── snippets
│   └── volumePartition
│       ├── volumePartition.png
│       └── volumePartition.py
├── testcases
│   ├── kelvin-helmholtz
│   │   ├── densityPlotter.py
│   │   └── generateIC.py
│   └── sedov
│       ├── initial_sedov.py
│       └── sedov_N61.h5
└── tools
    └── plotInteraction.py
```
### Directories

- `demonstrator` is dedicated to write a single CPU C++ program to implement the algorithm
- `media` holds resources for readme files
- `snippets` is for miscellaneous scripts to fool around 
- `testcases` includes subdirectories with standard test cases to validate the code

## Demonstrator 

The C++ program in the folder demonstrator can be build with the given `Makefile`. To install the necessary header-only libraries you can run the command `make install`.

The following prerequisites must be installed manually:

- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) for file-I/O
- [Boost property tree](https://www.boost.org/doc/libs/1_65_1/doc/html/property_tree.html) for configuration file parsing


## Visualization

### Volume partition

The fraction of volume $\psi_i\left(\vec{x}\right)$ assoziated with a particle $i$ is given by
$$\psi_i\left(\vec{x}\right) \equiv \frac{1}{\omega(\vec{x})}W\left(\vec{x}-\vec{x}_i, h\left(\vec{x}\right)\right)$$

where
$$\omega\left(\vec{x}\right) = \sum_j W\left(\vec{x}-\vec{x}_j, h\left(\vec{x}\right)\right)$$

and $h\left(\vec{x}\right)$ is some kernel size of the kernel function $W$.

The volume partition for a random particle distribution and periodic boundary conditions and a cubic spline kernel is shown below:

![volumePartition](media/volumePartition.png)

The plot can be generated in `snippets/volumePartition`.

## Test problems

### Kelvin-Helmholtz instability

A 2D test case setup as described in [McNally et. al. (2012)](https://arxiv.org/abs/1111.1764). 

<!--

### Sedov blast wave


### Keplerian disk
-->