import numpy as np
# minimal implementation with One square (four nodes)
"""
N4 —— N3
|     |
|     |
|     |
N1 —— N2

First triangle: globals N1, N2, N4
Second triangle: globals N2, N3, N4 -- n1, n2, n3 local

"""
n_triangles = 2
n_nodes= ((n_triangles*2)+1)**2 # (n_squares+ 1)**2 number of global nodes in the mesh



# NOTE test values
k= 1 # TODO input file
h = 1 # TODO input file -- h is the element thickness

class Node:
    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y
class Triangle: # element class
    def __init__(self, id, n1, n2, n3, k=1, h=1):
        self.id= id
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.k = k #TODO input file
        self.h = h #TODO input file
    def area(self):
        x1, y1 = self.n1.x, self.n1.y
        x2, y2 = self.n2.x, self.n2.y
        x3, y3 = self.n3.x, self.n3.y
        return 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
    def a(self, i, j):
        return i.x*j.y -j.x*i.y
    def b(self, i, j):
        """X coord diff between each node of the triangle"""
        return i.y - j.y
    def c(self, i, j):
        """Y coord diff between each node of the triangle"""
        return i.x - j.x
    def a_coeffs(self):
        """Area coefficients for the triangle"""
        return [self.a(self.n2, self.n3), self.a(self.n3, self.n1), self.a(self.n1, self.n2)]
    def b_coeffs(self): 
        """b coefficients for the triangle"""
        return [self.b(self.n2, self.n3), self.b(self.n3, self.n1), self.b(self.n1, self.n2)]
    def c_coeffs(self): 
        """c coefficients for the triangle"""
        return [self.c(self.n3, self.n2), self.c(self.n1, self.n3), self.c(self.n2, self.n1)]
    def local_stiffness_matrix(self): #CHECK
        # k is element dependent 
        A = self.area()
        k = self.k
        h = self.h
        b_coeffs = self.b_coeffs()
        c_coeffs = self.c_coeffs()
        for i in range(3):
            for j in range(3):
                k_local[i][j] = (b_coeffs[i] * b_coeffs[j] + c_coeffs[i] * c_coeffs[j]) * (k*h) / (4 * A)
        return k_local
class Mesh:
    def __init__(self, n_triangles):
        self.n_triangles = n_triangles  # total number of triangles in the mesh
        self.nodes = []
        self.triangles = []
        self.create_triangle_mesh()

    def create_triangle_mesh(self):
        n_rows = int((self.n_triangles / 2) ** 0.5)  # number of rows of triangles
        for i in range(n_rows + 1):  # number of nodes in each row
            for j in range(n_rows + 1):  # number of nodes in each column
                self.add_node(j, i)  # Create nodes directly with (x, y) coordinates

        for i in range(n_rows):  # Create triangles directly
            for j in range(n_rows):
                # Global node indices for the triangle
                n1_idx = i * (n_rows + 1) + j
                n2_idx = n1_idx + 1
                n3_idx = n1_idx + (n_rows + 1)
                n4_idx = n3_idx + 1

                # First triangle: N1, N2, N3
                self.triangles.append(Triangle(len(self.triangles), self.nodes[n1_idx], self.nodes[n2_idx], self.nodes[n3_idx]))

                # Second triangle: N2, N4, N3
                self.triangles.append(Triangle(len(self.triangles), self.nodes[n2_idx], self.nodes[n3_idx], self.nodes[n4_idx]))


            # DEBUG global node numbers for the square
                print(f"Square ({i}, {j}): N1={n1_idx} ({self.nodes[n1_idx].x}, {self.nodes[n1_idx].y}), "
                    f"N2={n2_idx} ({self.nodes[n2_idx].x}, {self.nodes[n2_idx].y}), "
                    f"N3={n3_idx} ({self.nodes[n3_idx].x}, {self.nodes[n3_idx].y}), "
                    f"N4={n4_idx} ({self.nodes[n4_idx].x}, {self.nodes[n4_idx].y})")

            # Plot the square and its triangles for visualization
                import matplotlib.pyplot as plt

                # Plot the current square
                square_x = [self.nodes[n1_idx].x, self.nodes[n2_idx].x, self.nodes[n4_idx].x, self.nodes[n3_idx].x, self.nodes[n1_idx].x]
                square_y = [self.nodes[n1_idx].y, self.nodes[n2_idx].y, self.nodes[n4_idx].y, self.nodes[n3_idx].y, self.nodes[n1_idx].y]
                plt.plot(square_x, square_y, 'b-', label='Square' if len(self.triangles) == 0 else "")

                # Plot the triangles within the square
                tri1 = self.triangles[-2]  # First triangle of the square
                tri2 = self.triangles[-1]  # Second triangle of the square

                tri1_x = [tri1.n1.x, tri1.n2.x, tri1.n3.x, tri1.n1.x]
                tri1_y = [tri1.n1.y, tri1.n2.y, tri1.n3.y, tri1.n1.y]
                plt.plot(tri1_x, tri1_y, 'r-', label='Triangle 1' if len(self.triangles) == 2 else "")

                tri2_x = [tri2.n1.x, tri2.n2.x, tri2.n3.x, tri2.n1.x]
                tri2_y = [tri2.n1.y, tri2.n2.y, tri2.n3.y, tri2.n1.y]
                plt.plot(tri2_x, tri2_y, 'g-', label='Triangle 2' if len(self.triangles) == 2 else "")
                plt.legend()

                # Add labels for nodes
                for idx, node in enumerate(self.nodes):
                    plt.text(node.x, node.y, f'N{idx}', fontsize=8, ha='right')

    def add_node(self, x, y):
        """
        Add a node to the mesh if it doesn't already exist.
        Return the existing or newly created node.
        """
        for node in self.nodes:
            if node.x == x and node.y == y:
                return node
        new_node = Node(x, y)
        self.nodes.append(new_node)
        return new_node

    def get_nodes(self):
        """
        Return all unique nodes in the mesh.
        """
        return self.nodes

    def get_triangles(self):
        """
        Return all triangles in the mesh.
        """
        return self.triangles

    def get_global_node(self, index):
        """
        Retrieve a global node by its index.
        """
        if 0 <= index < len(self.nodes):
            node = self.nodes[index]
            print(f"Node {index}: (x={node.x}, y={node.y})")
            return node
        else:
            raise IndexError("Global node index out of range.")

