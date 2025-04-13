from min_FEM import *
import numpy as np
import scipy.integrate

def compute_local_stiffness_matrix(k, element_id, mesh):
    # k: heat conductivity coefficient
    # tri_id: id of the triangle in the mesh
    # mesh: mesh object

    H_e = np.zeros((3, 3)) # local stiffness matrix

    # get element
    element = mesh.elements[element_id-1]
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

def compute_free_nodes(H, mesh):
    dirichlet_nodes = mesh.dirichlet_nodes
    # matrix without Dirichlet nodes (free nodes)
    H_free = np.delete(H, [node.id-1 for node in dirichlet_nodes], axis=0)
    H_free = np.delete(H_free, [node.id-1 for node in dirichlet_nodes], axis=1)

    return H_free

def compute_load_vector(H, mesh, q):
    # f = - integral(N_a * q)
    f = np.zeros(len(mesh.nodes)) # global load vector
    neumann_nodes = mesh.neumann_nodes
    elements = mesh.elements
    for element in elements:
        f_e = np.zeros(3) # local load vector
        a1, a2, a3 = element.a_coeffs()
        b1, b2, b3 = element.b_coeffs()
        c1, c2, c3 = element.c_coeffs()

        # basis functions
        A = element.area()
        N1 = lambda x, y: 1/(2*A) * (a1 + b1*x + c1*y)
        N2 = lambda x, y: 1/(2*A) * (a2 + b2*x + c2*y)
        N3 = lambda x, y: 1/(2*A) * (a3 + b3*x + c3*y)

        # find edges on the neumann boundary
        if element.n1 in neumann_nodes and element.n2 in neumann_nodes:
            # edge n1-n2
                
            # edge length
            l = np.sqrt((element.n1.x - element.n2.x)**2 + (element.n1.y - element.n2.y)**2) # length of the edge'

            # integral over N_a*q = l/2 * q
            f_e[0] = l/2 * q * N1(element.n1.x, element.n1.y) # N1 is the basis function for node n1
            f_e[1] = l/2 * q * N2(element.n2.x, element.n2.y) # N2 is the basis function for node n2
            f_e[2] = 0 # N3 is not on the edge            

        elif element.n2 in neumann_nodes and element.n3 in neumann_nodes:
            # edge n2-n3
            l = np.sqrt((element.n2.x - element.n3.x)**2 + (element.n2.y - element.n3.y)**2)
            f_e[0] = 0 # N1 is not on the edge
            f_e[1] = l/2 * q * N2(element.n2.x, element.n2.y) # N2 is the basis function for node n2
            f_e[2] = l/2 * q * N3(element.n3.x, element.n3.y) # N3 is the basis function for node n3
           
        elif element.n3 in neumann_nodes and element.n1 in neumann_nodes:
            # edge n3-n1
            l = np.sqrt((element.n3.x - element.n1.x)**2 + (element.n3.y - element.n1.y)**2)
            f_e[0] = l/2 * q * N1(element.n1.x, element.n1.y) # N1 is the basis function for node n1
            f_e[1] = 0 # N2 is not on the edge
            f_e[2] = l/2 * q * N3(element.n3.x, element.n3.y) # N3 is the basis function for node n3
            
        f[element.n1.id-1] += f_e[0]
        f[element.n2.id-1] += f_e[1]
        f[element.n3.id-1] += f_e[2]

    return -f



def compute_rhs(H, mesh, f, T_dirichlet):
    neumann_nodes = mesh.neumann_nodes
    dirichlet_nodes = mesh.dirichlet_nodes
    neumann_nodes_inside = mesh.neumann_nodes_inside
    free_nodes = np.concatenate([neumann_nodes_inside, neumann_nodes]) 

    
    # get vector P_11...100 for right hand side
    P_1 = f[[node.id -1 for node in free_nodes]]
    
    # get submatrix of H for Neumann nodes (H_1...10,11...100)
    H_neumann = np.delete(H, [node.id-1 for node in dirichlet_nodes], axis=0) # removed rows corresponding to Dirichlet nodes 
    H_neumann = np.delete(H_neumann, [node.id-1 for node in free_nodes], axis=1) # removed columns corresponding to free nodes
    
    T_0 = np.full(len(dirichlet_nodes), T_dirichlet)

    rhs = P_1 - np.matmul(H_neumann, T_0) # right hand side of the system of equations
    
    return rhs

def compute_reaction_forces(H, mesh, T):
    P = np.zeros(len(mesh.nodes)) # initialization of T_1...100
    P = np.matmul(H, T) 

    dirichlet_nodes = mesh.dirichlet_nodes
    return P[[node.id-1 for node in dirichlet_nodes]]# reaction forces at Dirichlet nodes
    