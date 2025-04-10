import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from min_FEM import *

class Mesh:
    def __init__(self, N):
        self.N = N
        self.nodes = []  # Initialize node list
        self.elements = []

        self.create_nodes(N)
        self.create_elements(N)

        
        self.neumann_nodes, self.dirichlet_nodes, self.neumann_nodes_inside = self.boundary_conditiones()
    
    
    def coordinates(self):
        L = 1.0  # Length of the domain
        # Create a grid of points
        x = np.linspace(0, L, int(math.sqrt(self.N)))
        y = np.linspace(0, L, int(math.sqrt(self.N)))
        X, Y = np.meshgrid(x, y)
        return X.flatten(), Y.flatten()
    
    def transform():
        return None
    
    def create_nodes(self, N):
        # self.node_list = []
        n = int(math.sqrt(N))
        # Create nodes with correct numbering
        x_coords, y_coords = self.coordinates()
        for i in range(N):
            x = x_coords[i]
            y = y_coords[i]
            self.nodes.append(Node(i+1, x, y))

    
    def get_node(self, index):
        return self.nodes[index - 1]  # Convert 1-based to 0-based indexing
    
    def create_elements(self, N):  # Fixed syntax and parameters
        n = int(math.sqrt(N))
        # Create elements with proper connectivity
        for j in range(n-1):
            for i in range(n-1):
                node1 = self.get_node(j * n + i + 1)
                node2 = self.get_node(j * n + i + 2)
                node3 = self.get_node((j + 1) * n + i + 1) 
                node4 = self.get_node((j + 1) * n + i + 2)
                # Create two triangular elements
                self.elements.append(Triangle(j*n+i+1, node1, node2, node3))  # Lower triangle
                self.elements.append(Triangle(j*n+i+1, node2, node4, node3))  # Upper triangle

    
    def boundary_conditiones(self):
        # Define Dirichlet and Neumann boundary conditions
        dirichlet_nodes = (self.node_list[i] for i in range(1, 11)) 
        # Neumann BCs (N91 to N100)
        neumann_nodes= (self.node_list[i] for i in range(91, 101))
        neumann_nodes_inside= (self.node_list[i] for i in range(11, 91))
        return dirichlet_nodes, neumann_nodes, neumann_nodes_inside
    
    def plot_mesh(self):
        # Get coordinates
        x, y = self.coordinates()

        # Create triangulation
        triang = Triangulation(x, y, self.elements)

        # Plot the mesh
        plt.figure(figsize=(8, 8))
        plt.triplot(triang, color='gray', alpha=0.5)
        plt.scatter(x, y, color='red')
        plt.title('Mesh Plot')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.show()

           



