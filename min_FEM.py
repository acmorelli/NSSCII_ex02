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
n_square= 1
n_nodes= (n_square**2 +1)**2
n_triangles = n_square * 2

# NOTE test values
k= 1 # TODO input file
h = 1 # TODO input file -- h is the element thickness

class Node:
    def __init__(self, x, y):
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
    def b(self, i, j):
        """X coord diff between each node of the triangle"""
        return i.x - j.x
    def c(self, i, j):
        """Y coord diff between each node of the triangle"""
        return i.y - j.y
    def b_coeffs(self): 
        """b coefficients for the triangle"""
        return [self.b(self.n1, self.n2), self.b(self.n2, self.n3), self.b(self.n3, self.n1)]
    def c_coeffs(self): 
        """c coefficients for the triangle"""
        return [self.c(self.n1, self.n2), self.c(self.n2, self.n3), self.c(self.n3, self.n1)]
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
class Square:
    def __init__(self, id, n1, n2, n3, n4):
        self.id = id
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.n4 = n4
    def triangles(self):
        return [Triangle(self.id, self.n1, self.n2, self.n4), Triangle(self.id, self.n2, self.n3, self.n4)]
    
class Mesh: #TODO
    def __init__(self, n_square):
        self.n_square = n_square
        self.squares = create_square_mesh(n_square) #TODO
        self.triangles = []
        for square in self.squares:
            self.triangles.extend(square.triangles())

    def create_square_mesh(self):
        """
        square mesh n_square x n_square squares
        """
        squares = []
        for i in range(self.n_square):
            for j in range(self.n_square):
                n1 = Node(i, j)
                n2 = Node(i+1, j)
                n3 = Node(i+1, j+1)
                n4 = Node(i, j+1)
                self.square = Square(len(self.squares), n1, n2, n3, n4)
                self.squares.append(self.square)
        return self.squares

    def get_nodes(self):
        nodes = []
        for square in self.squares:
            nodes.extend([square.n1, square.n2, square.n3, square.n4])
        return nodes
    def get_triangles(self):
        triangles = []
        for square in self.squares:
            triangles.extend(square.triangles())


    def mapping_global(self):
        """
        Map the global node numbers to the local node numbers
        """
        return [self.n1, self.n2, self.n3, self.n4]
    def mapping_local(self):
        """
        Map the local node numbers to the global node numbers

        """
        return [self.n1, self.n2, self.n3, self.n4]

"""
# Workflow:

Create mesh of squares and triangles
Compute area of each triangle
Compute b and c coefficients (based on coordinates of nodes) of each triangle
Compute local stiffness matrix of each triangle


next
mesh class -- create mesh properly
connectivity matrix -- local to global
create global stiffness matrix

force vector

"""



