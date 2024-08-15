# LevelSet-VOF
VOF+LevelSet method implementation in OpenFoam - Currently Supported for OpenFOAM V6.0
Contact : amirhossein.taran@ucdconnect.ie
______________________________________________
## Introduction
The level-set method was developed in 1979 by Alain Dervieux, and subsequently popularized by Stanley Osher and James Sethian. To run a Math Model in the interface of two different fluids we need to soften the interactions between the fluids. Therefore we need to apply a specific function: Compact Level Set Method.
As a “spin off”, the CompactLSM is a complement of the LSM, that helps solving LSM equations. It can be used in numerical simulation of flow, for example, if we are working with discretization of the interface water-air, compacts at sixth order, ensures the accurate and fast calculation of the interface equations (Monteiro 2018).
The LSM uses a distance function to locate different fluids. A distance function is that whose value represents the smallest distance from the point where it is being analyzed to the interface. This distance function is identified by isolines (2D) or isosurfaces (3D), showing that  the negative values refer to one of the fluids, positive values refer to the other and the zero value corresponds to the position of the interface.
Since the specific mass and viscosity are discontinuous at the interface, both excess diffusion problem (interface widening) and numerical oscillations are expected if there is no adequate treatment of the fluid near the interface. To minimize these problems, the Level Set method uses a smooth, cell-related Heaviside function that explicitly defines the interface position.

![grab-landing-page](https://github.com/cloner0110/LevelSet-VOF/blob/main/ezgif-1-c13e0ed189.gif)



