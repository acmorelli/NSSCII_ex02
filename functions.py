from min_FEM import *
import numpy as np

def compute_local_stiffness_matrix(k, tri_id, mesh):
    # k: heat conductivity coefficient
    # tri_id: id of the triangle in the mesh
    # mesh: mesh object

    H_e = np.zeros((3, 3)) # local stiffness matrix

    # get element
    element = mesh.triangles[tri_id]
    # coefficients
    b = np.array(element.b_coeffs())
    c = np.array(element.c_coeffs())
    # Area of the triangle
    A = element.area()
    # see Zienkewicz, page 120 and 125
    for i in range(3):
        for j in range(3):
            H_e[i][j] = (k / (2*A)) * (b[i]*b[j] + c[i]*c[j]) 

    return H_e