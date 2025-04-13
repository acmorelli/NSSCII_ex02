from min_FEM import *
import numpy as np

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
    neumann_nodes = mesh.neumann_nodes
    elements = mesh.elements
    for element in elements:
        pass
    pass



def compute_rhs(H, mesh, T_dirichlet):
    neumann_nodes = mesh.neumann_nodes
    dirichlet_nodes = mesh.dirichlet_nodes
    neumann_nodes_inside = mesh.neumann_nodes_inside
    
    # get vector P_11...100 for right hand side
    P_1 = np.zeros(len(mesh.nodes) - len(dirichlet_nodes))
    
    # get submatrix of H for Neumann nodes (H_1...10,11...100)
    H_neumann = np.delete(H, [node.id-1 for node in dirichlet_nodes], axis=0) # removed rows corresponding to Dirichlet nodes 
    free_nodes = np.concatenate([neumann_nodes_inside, neumann_nodes]) 
    H_neumann = np.delete(H_neumann, [node.id-1 for node in free_nodes], axis=1) # removed columns corresponding to free nodes
    
    T_0 = np.full(len(dirichlet_nodes), T_dirichlet)

    rhs = P_1 - np.dot(H_neumann, T_0) # right hand side of the system of equations
    
    return rhs

def compute_reaction_forces(H, mesh, T_free):
    P = np.zeros(len(mesh.nodes)) # initialization of P_1...100
    T = np.zeros(len(mesh.nodes)) # initialization of T_1...100
    
    # intialization of Temperature vector with solution T_free
    neumann_nodes = mesh.neumann_nodes
    neumann_nodes_inside = mesh.neumann_nodes_inside
    free_nodes = np.concatenate([neumann_nodes_inside, neumann_nodes])
    for i, node in enumerate(free_nodes):
        T[node.id-1] = T_free[i] 
    
    P = np.dot(H, T) 

    
    dirichlet_nodes = mesh.dirichlet_nodes
    return P[[node.id for node in dirichlet_nodes]]# reaction forces at Dirichlet nodes
    