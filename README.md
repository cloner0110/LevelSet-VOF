# LevelSet-VOF
VOF+LevelSet method implementation in OpenFoam - V-6.0
______________________________________________
The level-set method was developed in 1979 by Alain Dervieux, and subsequently popularized by Stanley Osher and James Sethian. To run a Math Model in the interface of two different fluids we need to soften the interactions between the fluids. Therefore we need to apply a specific function: Compact Level Set Method.
As a “spin off”, the CompactLSM is a complement of the LSM, that helps solving LSM equations. It can be used in numerical simulation of flow, for example, if we are working with discretization of the interface water-air, compacts at sixth order, ensures the accurate and fast calculation of the interface equations (Monteiro 2018).
The LSM uses a distance function to locate different fluids. A distance function is that whose value represents the smallest distance from the point where it is being analyzed to the interface. This distance function is identified by isolines (2D) or isosurfaces (3D), showing that  the negative values refer to one of the fluids, positive values refer to the other and the zero value corresponds to the position of the interface.
Since the specific mass and viscosity are discontinuous at the interface, both excess diffusion problem (interface widening) and numerical oscillations are expected if there is no adequate treatment of the fluid near the interface. To minimize these problems, the Level Set method uses a smooth, cell-related Heaviside function that explicitly defines the interface position (![equation](http://www.sciweavers.org/upload/Tex2Img_1652109100/render.png)).

The initial value of Level set function is given by:

![equation](http://www.sciweavers.org/upload/Tex2Img_1652109424/render.png)

by using this initial value, and solving Eikonal Equation, the equations become :

![equation](http://www.sciweavers.org/upload/Tex2Img_1652109874/render.png)

in this equation, ![equation](http://www.sciweavers.org/upload/Tex2Img_1652110030/render.png) is an artificial time step and needs to small enough to minimze changes in Levelset function, this value commonly considered as:

![equation](http://www.sciweavers.org/upload/Tex2Img_1652110131/render.png)

finally, the LevelSet equation obtained as follows:

![equation](http://www.sciweavers.org/upload/Tex2Img_1652110278/render.png)

by obtaining LevelSet function from previous equation, the smoothed Heaviside function can be introduced as follows:

![equation](http://www.sciweavers.org/upload/Tex2Img_1652110555/render.png)

at last, the surface tension force is obtained using mentioned Heaviside function:

![equation](http://www.sciweavers.org/upload/Tex2Img_1652110686/render.png)



