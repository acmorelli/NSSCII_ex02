import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from min_FEM import *

class Mesh:
    def __init__(self, N, Variation, L):
        self.N = N
        self.nodes = []  # Initialize node list
        self.elements = []
        self.Variation = Variation
        self.create_nodes(N)
        self.create_elements(N)
        self.L=L

        
        self.neumann_nodes, self.dirichlet_nodes, self.neumann_nodes_inside = self.boundary_conditiones()
    
    
    def coordinates(self):
        # Create a grid of points
        x = np.linspace(0, self.L, int(math.sqrt(self.N)))
        y = np.linspace(0, self.L, int(math.sqrt(self.N)))
        X, Y = np.meshgrid(x, y)
        return X.flatten(), Y.flatten()
 
        
    
    def create_nodes(self, N):
        self.Variation
        # self.node_list = []
        x_bias=[] 
        n = int(math.sqrt(N))
        # Create nodes with correct numbering
        x_coords, y_coords = self.transform(self.Variation)
        for i in range(N):
            x = x_coords[i]
            y = y_coords[i]
            print("x", x)
            print("y", y)
            if  self.Variation == 'V1':
                if i == 9:
                    # Apply transformation for V1
                    x1 = self.L/2
                    y1 = 0
                    self.nodes.append(Node(i+1, x1, y1))
                else:
                    self.nodes.append(Node(i+1, x, y))
                # Apply transformation for V1   
            elif self.Variation == 'V2':
                # Define the bias factor B
                B = (1 / (2 * self.L)) * (self.L - y_coords[i])
                # Apply the bias to the x-coordinates
                x_bias[i] = x_coords[i] * ((B / self.L) * x_coords[i] - B + 1)
                self.nodes.append(Node(i+1, x_bias[i], y_coords))

            elif self.Variation == 'V3':
                r_min=0.5*self.L
                r_max=0.*5*self.L
                theta_min,theta_max=0,np.pi/4
                


                self.nodes.append(Node(i+1,r,theta))
            
           
            else:
                #Regular Mesh   
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

           



