# import local files here
from functions import compute_local_stiffness_matrix
from mesh import Mesh
import numpy as np

testing = True

def main():
    
    k = 0.1  # Example conductivity coefficient
    N = 100 # Number of nodes in the mesh (4x4 grid)

    assert np.sqrt(N) % 1 == 0, "N must be a perfect square"

    if testing:
        # Create 2 by 2 mesh
        N = 4
        mesh = Mesh(N)
        for node in mesh.nodes:
            print(f"Node: {node.id}, Coordinates: ({node.x}, {node.y})")
        for element in mesh.elements:
            print(f"Element: {element.id}, Nodes: ({element.n1.id}, {element.n2.id}, {element.n3.id})")

    else:
        # create real mesh
        pass
    
    
    
    H = np.zeros((N, N))  # Global stiffness matrix initialization
    for element_id in range(1, len(mesh.elements) + 1):

        # Compute local stiffness matrix for each element
        H_e = compute_local_stiffness_matrix(k, element_id, mesh)
        
        # Assemble global stiffness matrix
        element_nodes_ids = [mesh.elements[element_id-1].n1.id, mesh.elements[element_id-1].n2.id, mesh.elements[element_id-1].n3.id]
        print("Element id", element_id, "element_nodes", element_nodes_ids)
        for row_e in range(H_e.shape[0]):
            global_row_idx = element_nodes_ids[row_e] - 1
            for col_e in range(H_e.shape[1]):
                global_col_idx = element_nodes_ids[col_e] - 1
                print("global_row_idx", global_row_idx, "global_col_idx", global_col_idx)
                H[global_row_idx][global_col_idx] += H_e[row_e][col_e]
                

    
    



if __name__ == "__main__":
    main()