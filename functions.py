from min_FEM import *
import numpy as np

def compute_a(element: Triangle):
    # see definition of a in Zienkewicz, page 120
    a = np.zeros(3)
    n1 = element.n1
    n2 = element.n2
    n3 = element.n3
    a[0] = n2.x*n3.y - n3.x*n2.y
    a[1] = n3.x*n1.y - n1.x*n3.y
    a[2] = n1.x*n2.y - n2.x*n1.y
    return a

def compute_b(element: Triangle):
    # see definition of b in Zienkewicz, page 120
    b = np.zeros(3)
    n1 = element.n1
    n2 = element.n2
    n3 = element.n3
    b[0] = n2.y - n3.y
    b[1] = n3.y - n1.y
    b[2] = n1.y - n2.y
    return b

def compute_c(element: Triangle):
    # see definition of c in Zienkewicz, page 120
    c = np.zeros(3)
    n1 = element.n1
    n2 = element.n2
    n3 = element.n3
    c[0] = n3.x - n2.x
    c[1] = n1.x - n3.x
    c[2] = n2.x - n1.x
    return c

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
            H_e[i][j] = (k / (2*A)) *(b[i]*b[j] + c[i]*c[j]) 

    return H_e


def compute_local_force_vector(mesh):
    pass

def compute_connectivity_matrix(mesh):
    pass