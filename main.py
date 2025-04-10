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
    



if __name__ == "__main__":
    main()