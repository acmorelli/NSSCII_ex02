# import local files here
from functions import *
from mesh import Mesh
import numpy as np

testing = True

def main():
    
    k = 1  # Example conductivity coefficient
    N = 100 # Number of nodes in the mesh (4x4 grid)
    T_dirichlet = 1  # Dirichlet boundary condition value
    q = 1  # Flux across the Neumann boundary

    assert np.sqrt(N) % 1 == 0, "N must be a perfect square"

    if testing:
        # Create 2 by 2 mesh
        N = 4
        L = 1
        mesh = Mesh(N,'V3', L)
        mesh.plot_mesh()
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
        H_e = compute_local_stiffness_matrix(k, element_id, mesh)
        
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
    rhs = compute_rhs(H, mesh, T_dirichlet)
    print("rhs:", rhs)
    # solve the system of equations
    T_free = np.linalg.solve(H_free, rhs)
    print("T_free:", T_free)
    reaction_forces = compute_reaction_forces(H, mesh, T_free)
    print("Reaction forces:", reaction_forces)
    # post processing
        
                

    
    



if __name__ == "__main__":
    main()