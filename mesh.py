import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

class Mesh:
    #TO DO 
    #1-->100 points  YES
    #2-->N1(0,0), N10(L,0),N100,(L,L) YES
    #3-->Element connectivity scheme  where every element En is assigned such as E1(N1,N2,N11)
    #Counterclockwise nodal assignment has to be used YES

    #-->Boundary conditions: N1- N10 Dirichlet BCs
    #-->N91 – N100: Neumann BCs 
    #-->o,w,Neumann BCs P=0
    
    #IS IT COMPLETED UNTIL HERE? :  NO

    #6 Variations of mesh also have to be solved

    # After Mesh creation one have to divide L(length on both x and y axis) and hz(tickness at z axis) 
    # to meshing for group specific work

    def __init__(self, N):
        self.N = N
        self.node_list = []
        self.elements=[]  # Initialize node list
    
    def coordinates(self):
        L = 1.0  # Length of the domain
        # Create a grid of points
        x = np.linspace(0, L, int(math.sqrt(self.N)))
        y = np.linspace(0, L, int(math.sqrt(self.N)))
        X, Y = np.meshgrid(x, y)
        return X.flatten(), Y.flatten()
    
    def transforrm():
        return None
    
    def nodes(self, N):
        self.node_list = []
        n = int(math.sqrt(N))
        # Create nodes with correct numbering
        for j in range(n):
            for i in range(n):
                node_num = j * n + i + 1  # 1-based numbering
                self.node_list.append(node_num)
        print(self.node_list)
        return self.node_list 
    
    def get_node(self, index):
        return self.node_list[index - 1]  # Convert 1-based to 0-based indexing
    
    def connectivity(self, N):  # Fixed syntax and parameters
        n = int(math.sqrt(N))
        # Create elements with proper connectivity
        for j in range(n-1):
            for i in range(n-1):
                node1 = self.get_node(j * n + i + 1)
                node2 = self.get_node(j * n + i + 2)
                node3 = self.get_node((j + 1) * n + i + 1) 
                node4 = self.get_node((j + 1) * n + i + 2)
                # Create two triangular elements
                self.elements.append([node1, node2, node3])  # Lower triangle
                self.elements.append([node2, node4, node3])  # Upper triangle
        print(self.elements)
        return self.elements
    
    def boundary_conditioned(self):
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

           



