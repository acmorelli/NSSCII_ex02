# import local files here
from functions import compute_local_stiffness_matrix
from mesh import Mesh

testing = True

def main():
    if testing:
        # Create 2 by 2 mesh
        mesh = Mesh(4)
        for node in mesh.nodes:
            print(f"Node: {node.id}, Coordinates: ({node.x}, {node.y})")
        for element in mesh.elements:
            print(f"Element: {element.id}, Nodes: ({element.n1.id}, {element.n2.id}, {element.n3.id})")

    else:
        # create real mesh
        pass
    
    # Compute local stiffness matrix for each element
    k = 0.1  # Example conductivity coefficient
    for element_id in range(1, len(mesh.elements) + 1):
        H_e = compute_local_stiffness_matrix(k, element_id, mesh)
        # continue code here



if __name__ == "__main__":
    main()