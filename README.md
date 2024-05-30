# MPhil Project

This repository contains all the code files for the MPhil Project titled "CFD Simulation in Geometries with Multiple Rotating Frames of Reference" by Alyssa Gomes.

The code generates a computational mesh using input data about points, faces, cells, boundaries, face zones, and cell zones. It is capable of discretizing terms like convection, diffusion (Laplacian), and temporal derivatives, thus being able to solve the incompressible Navier-Stokes equations, SRF and MRF models. 

The following are the external libraries and software packages used in this project:
- OpenFOAM: to generate mesh input files
- libconfig++: used for the settings files
- Eigen: used to store and solve the matrices used in the ScalarField and VectorField classes
- VTK9: used to output VTK XML files for post-processing using applications like ParaView or VisIt
- ParaView: To visualize the VTK output files
- CMake
- Makefile


To compile and build a specific project in Linux, run:
```bash
$ make clean
$ make build
$ cd build
$ make main_SRF
$ cd .. 
$ ./bin/main_SRF SRF.cfg
```

To visualize the results run either of the following:
```bash
$ paraview
$ visit
```

To clear the VTK output folder of a test case, run:
```bash
$ make clean-VTK test_case=mixerSRF
```
