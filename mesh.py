import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from min_FEM import *

class Mesh:
    def __init__(self, N, Variation, L,k,y_neumann, y_dirichlet, hz):
        self.N = N
        self.nodes = []  # Initialize node list
        self.elements = []
        self.Variation = Variation
        self.L = L #length
        self.hz = hz # thickness
        self.create_nodes(N, Variation, L)  # Pass all required parameters
        self.create_elements(N, k, hz)
        self.neumann_nodes, self.dirichlet_nodes, self.neumann_nodes_inside = self.boundary_conditiones(y_neumann, y_dirichlet)
    
    def coordinates(self):
        # Create a grid of points
        x = np.linspace(0, self.L, int(math.sqrt(self.N)))
        y = np.linspace(0, self.L, int(math.sqrt(self.N)))
        X, Y = np.meshgrid(x, y)
        return X.flatten(), Y.flatten()
 
        
    
    def create_nodes(self, N, Variation, L):  # Match parameters with actual usage
        x_bias = []  # Initialize x_bias list
        n = int(math.sqrt(N))

        #Creating grid of points for V3 variation
        r1 = np.linspace(L, 2*L, int(math.sqrt(self.N)))
        theta1 = np.linspace(np.pi,5*np.pi/4 , int(math.sqrt(self.N)))
        R, Theta = np.meshgrid(r1, theta1)
        R = R.flatten()
        Theta = Theta.flatten()

        # Create nodes with correct numbering
        x_coords, y_coords = self.coordinates()
        for i in range(N):
            r=R[i]
            theta=Theta[i]
            x = x_coords[i]
            y = y_coords[i]
            # print("x", x)
            # print("y", y)
            if Variation == 'V1':
                # Apply transformation for V1   
                x =  (1+y) * x/2
                self.nodes.append(Node(i+1, x, y))
                
            elif Variation == 'V2':
                # Define the bias factor B
                B = (1 / (2 * L)) * (L - y)
                # Apply the bias to the x-coordinates
                x_bias.append(x * ((B / L) * x - B + 1))
                self.nodes.append(Node(i+1, x_bias[i], y))

            elif Variation == 'V3':
                # Apply transformation for V3
                X_p = r * np.cos(theta) + 2*L
                Y_p = r * np.sin(theta)

                self.nodes.append(Node(i+1,X_p,Y_p))
               
            else:
                #Regular Mesh -- squared easy mesh???
                self.nodes.append(Node(i+1, x, y))

    
    def get_node(self, index):
        return self.nodes[index - 1]  # Convert 1-based to 0-based indexing
    
    def create_elements(self, N, k, hz):  # Fixed syntax and parameters
        n = int(math.sqrt(N))
        # Create elements with proper connectivity
        for j in range(n-1):
            for i in range(n-1):
                node1 = self.get_node(j * n + i + 1)
                node2 = self.get_node(j * n + i + 2)
                node3 = self.get_node((j + 1) * n + i + 1) 
                node4 = self.get_node((j + 1) * n + i + 2)
                
                # Create two triangular elements
                self.elements.append(Triangle((j*(n-1)+i)*2+1, node1, node2, node3, k, hz))  # Lower triangle
                self.elements.append(Triangle((j*(n-1)+i)*2+2, node2, node4, node3, k, hz))  # Upper triangle
                #print("Node1: ", node1.x, node1.y)
        # print("elements:",self.elements[0].n1.x,self.elements[0].n1.y)
        # print("elements:",self.elements[0].n2.x,self.elements[0].n2.y)
        # print("elements:",self.elements[0].n3.x,self.elements[0].n3.y)

    
    def boundary_conditiones(self, y_neumann, y_dirichlet):
        # Define Dirichlet and Neumann boundary conditions
        dirichlet_nodes = [node for node in self.nodes if node.y == y_dirichlet]
        neumann_nodes = [node for node in self.nodes if node.y == y_neumann]  
        neumann_nodes_inside = [node for node in self.nodes if node.id not in ({n.id for n in dirichlet_nodes} | {n.id for n in neumann_nodes})]        
        #assert len(neumann_nodes_inside) == len(self.nodes) - len(dirichlet_nodes) - len(neumann_nodes), "Mismatch in BCs nodes"
        return neumann_nodes, dirichlet_nodes, neumann_nodes_inside
    
    def plot_mesh(self):
        # Extract all node coordinates
        x_coords = np.array([node.x for node in self.nodes])
        y_coords = np.array([node.y for node in self.nodes])
        
        # Create triangulation connectivity list
        triangles = np.array([[self.nodes.index(elem.n1),
                             self.nodes.index(elem.n2),
                             self.nodes.index(elem.n3)] 
                            for elem in self.elements])
        
        # Create triangulation
        triang = Triangulation(x_coords, y_coords, triangles)

        # Plot the mesh
        plt.figure(figsize=(8, 8))
        plt.triplot(triang, color='gray', alpha=0.5)
        plt.scatter(x_coords, y_coords, color='red', s=20)
        plt.title('Mesh Plot')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.axis('equal')
        plt.show()
           


