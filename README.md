# Content of the files
## Overview
Our code is comprised of 4 files: main.py, min_FEM.py, mesh.py and function.py. Each file and its usage is explained below.
## main.py
This file is the main code that is to be executed. It contains the workflow of the FEM implementation. It defines the necessary parameters in the first lines of the main function. The parameters are the group specific input data as well as the variation ('V0', 'V1', 'V2', 'V3', 'V4a', 'V4b').
The las few lines of the code call the plotting functions, they can be commented out at will.

## min_FEM.py and mesh.py
The code has a object oriented implementation. These to file define the necessary classes, which are
### in min_FEM.py:
Node:
Defines a node with coordinates and node id.
Triangle:
Defines a triangle with with some useful functions and attributes like its id, the nodes its made of, geometric information like its area. Also FEM results like the flux will be stored in the Triangle objects.

### in mesh.py
Mesh:
This class defines the mesh. It contains all its Nodes and Elements as well as the boundary conditions.

## function.py
In order to keep the main.py clean and understandable, we created this file to contain all the complicated function. Particularly, they contain the functions that compute the stiffness matrix, force vector / right hand side as well as functions for post processing like the computation of the flux or the function for plotting the results.

# How to use
To use our code just execute the main.py. You can adjust the input parameters at the beginning of the main function and you can plot the desired output at the end of the main function by commenting out/in the respective functions. The parameter variables as well as the plotting function are clearly named.