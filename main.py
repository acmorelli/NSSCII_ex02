""" observations version 1 - squared mesh:
• The temperature gradient must be constant in the entire model and equal to the overall
gradient ΔT/Δy.
• The flux must be constant in the entire model and must be equal to the applied flux. The
latter must be equal to the sum of the applied nodal “forces” divided by the corresponding
area.
• The temperature gradient and the flux, together with eqn. (2) must yield the input for
the conductivity k.

SI-units to be used:
    + T in K
    + L in m
    + k in W/(mK)
    + q in W/m^2 - ad Neumann
    + P in W - ad nodal forces

# Length in x- and y-direction.
L = 0.1

# Thickness (z-direction).
hz = 0.01

# Thermal conductivity (k=k_xx=k_yy, k_xy = 0.).
k = 373.

# Factor c for modifying thermal conductivity k for
# elements in elements_to_be_modified.
c = 10.

# Elements to be modified.
elements_to_be_modified = [
                          62-70,
                          83-86,
                          98-103,
                          111-122
                          ]

# Boundary conditions.
q(y=0) = 1000000.
T(y=L) = 313.
"""
# import local files here
from functions import *
from mesh import Mesh
import numpy as np
import matplotlib.pyplot as plt
from print_HTP_2025 import print_HTP


testing = True

def main():
    
    k = 373  #  W/mKn
    # convert to W/K
    # k = k * 1e3 # W/K

    L = 0.1  # m (length of squared domain (V0))
    N = 100 # Number of nodes in x=0 (V0)
    q_neumann = 1000000  # w/ m2 (Flux across the Neumann boundary)
    y_neumann = 0.0 # y coordinate of the Neumann boundary
    T_dirichlet = 313.0 # K (Dirichlet bc)
    y_dirichlet = L # y coordinate of the Dirichlet boundary
    hz = 0.01 # m (thickness in z-direction)
    Variation=  'V0' # also change the mesh instantiation !!!
    
    """
    # # # debug
    L=1
    k=1
    q_neumann=1
    y_neumann= L
    T_dirichlet = 1
    y_dirichlet = 0.0
    """

    # # # # variation 4 # # # #
    id1=np.arange(62,71)
    id2=np.arange(83,87)
    id3=np.arange(98,104)
    id4=np.arange(111,123)
    id_c= np.concatenate((id1,id2,id3,id4))
    ce = 10.0 # factor for modifying thermal conductivity k 
    

    #assert np.sqrt(N) % 1 == 0, "N must be a perfect square"

    if testing:
        mesh = Mesh(N,Variation, L, k, y_neumann, y_dirichlet, hz) 

        for node in mesh.nodes:
            print(f"Node: {node.id}, Coordinates: ({node.x}, {node.y})")
        for element in mesh.elements:
            print(f"Element: {element.id}, Nodes: ({element.n1.id}, {element.n2.id}, {element.n3.id})")
        print("Dirichlet nodes:")
        for node in mesh.dirichlet_nodes:
            print(f"Node: {node.id}, Coordinates: ({node.x}, {node.y})")
        print("Neumann nodes:")
        for node in mesh.neumann_nodes:
            print(f"Node: {node.id}, Coordinates: ({node.x}, {node.y})")
        print("Neumann nodes inside:")
        for node in mesh.neumann_nodes_inside:
            print(f"Node: {node.id}, Coordinates: ({node.x}, {node.y})")

        print("\n")

    else:
        # create real mesh
        pass
    
    
    
    H = np.zeros((N, N))  # Global stiffness matrix initialization
    for element_id in range(1, len(mesh.elements) + 1):

        # Compute local stiffness matrix for each element
        H_e = compute_local_stiffness_matrix(k, element_id, mesh, hz, Variation, id_c, ce)
        
        # Assemble global stiffness matrix
        element_nodes_ids = [mesh.elements[element_id-1].n1.id, mesh.elements[element_id-1].n2.id, mesh.elements[element_id-1].n3.id]
        for row_e in range(H_e.shape[0]):
            global_row_idx = element_nodes_ids[row_e] - 1
            for col_e in range(H_e.shape[1]):
                global_col_idx = element_nodes_ids[col_e] - 1
                H[global_row_idx][global_col_idx] += H_e[row_e][col_e]

    # take into account the Dirichlet boundary conditions
    print("H:", H)
    H_free = compute_free_nodes(H, mesh)
    print("H_free:", H_free)
    f = compute_load_vector(H, mesh, q_neumann, hz)
    # # # # # 
    # # FIX FROM GPT
    # # Extract coupling matrix between free and Dirichlet nodes
    # dirichlet_ids = [node.id - 1 for node in mesh.dirichlet_nodes]
    # T_d = np.full(len(dirichlet_ids), T_dirichlet)

    # # Free nodes (complement of Dirichlet)
    # all_ids = np.arange(len(mesh.nodes))
    # free_ids = np.setdiff1d(all_ids, dirichlet_ids)

    # # Partition matrix
    # H_ff = H[np.ix_(free_ids, free_ids)]
    # H_fd = H[np.ix_(free_ids, dirichlet_ids)]

    # # Adjust RHS
    # f_free = f[free_ids] - H_fd @ T_d

    # # Solve
    # T_free = np.linalg.solve(H_ff, f_free)

    # # Reconstruct full solution
    # T = np.zeros(len(mesh.nodes))
    # T[dirichlet_ids] = T_dirichlet
    # T[free_ids] = T_free

    # print('T from gpt', T)

    
    print("f:", f)
    rhs = compute_rhs(H, mesh, f, T_dirichlet)
    print("rhs:", rhs)
    # solve the system of equations
    T_free = np.linalg.solve(H_free, rhs)
    print("T_free:", T_free)
    # set up complete solution vector
    T = np.zeros(len(mesh.nodes)) 
    # free_nodes = np.concatenate([mesh.neumann_nodes_inside, mesh.neumann_nodes])
    free_nodes = np.array([node for node in mesh.nodes if node not in mesh.dirichlet_nodes])
    for i, node in enumerate(free_nodes):
        T[node.id - 1] = T_free[i]
        print('node', node, 't_free[i]',T_free[i])
    for node in mesh.dirichlet_nodes:
        T[node.id - 1] = T_dirichlet
    print("T:", T)
    

    reaction_forces = compute_reaction_forces(H, mesh, T)
    print("Reaction forces:", reaction_forces)

    # add reaction forces to force vector
    counter = 0
    for node in mesh.nodes:
        if node in mesh.dirichlet_nodes:
            f[node.id - 1] += reaction_forces[counter]
            counter += 1
    print("Updated force vector:", f)

    # post processing
    compute_temperature_gradient(mesh, T)
    compute_heat_flux(mesh)


    # plotting
    plot_temperature_field(mesh, T, Variation)
    # plot_temperature_gradient(mesh, Variation)
    # plot_heat_flux(mesh, Variation)
    
    # save results to file
    filename = "output_" + Variation + ".txt"
    print_HTP(H, np.reshape(T, (-1,1)), np.reshape(f, (-1,1)), filename)
    print(f"Results saved to {filename}")

if __name__ == "__main__":
    main()