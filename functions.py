from min_FEM import *
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri

def compute_local_stiffness_matrix(k,element_id, mesh, Variation=None):

    # k: heat conductivity coefficient
    # tri_id: id of the triangle in the mesh
    # mesh: mesh object
    id1=np.arange(62,71)
    id2=np.arange(83,87)
    id3=np.arange(98,104)
    id4=np.arange(111,123)
    id_c=np.concatenate((id1,id2,id3,id4))
    ce=10
    H_e = np.zeros((3, 3)) # local stiffness matrix


    # get element
    element = mesh.elements[element_id-1]
    # coefficients
    b = np.array(element.b_coeffs())
    c = np.array(element.c_coeffs())
    # Area of the triangle
    A = element.area()
    if Variation == 'V4a':
        if element.id in id_c:
            k=ce*k
        else:
            pass
    elif Variation == 'V4b':
        if element.id in id_c:
            k=k/ce
        else:
            pass
    else:
        pass
    # see Zienkewicz, page 120 and 125
    for i in range(3):
        for j in range(3):
            H_e[i][j] = (k / (2*A)) * (b[i]*b[j] + c[i]*c[j]) 

    return H_e

def compute_free_nodes(H, mesh):
    dirichlet_nodes = mesh.dirichlet_nodes
    # matrix without Dirichlet nodes (free nodes)
    H_free = np.delete(H, [node.id-1 for node in dirichlet_nodes], axis=0)
    H_free = np.delete(H_free, [node.id-1 for node in dirichlet_nodes], axis=1)

    return H_free

def compute_load_vector(H, mesh, q):
    # f = - integral(N_a * q)
    f = np.zeros(len(mesh.nodes)) # global load vector
    neumann_nodes = mesh.neumann_nodes
    elements = mesh.elements
    for element in elements:
        f_e = np.zeros(3) # local load vector
        a1, a2, a3 = element.a_coeffs()
        b1, b2, b3 = element.b_coeffs()
        c1, c2, c3 = element.c_coeffs()

        # basis functions
        A = element.area()
        N1 = lambda x, y: 1/(2*A) * (a1 + b1*x + c1*y)
        N2 = lambda x, y: 1/(2*A) * (a2 + b2*x + c2*y)
        N3 = lambda x, y: 1/(2*A) * (a3 + b3*x + c3*y)

        # find edges on the neumann boundary
        if element.n1 in neumann_nodes and element.n2 in neumann_nodes:
            # edge n1-n2
                
            # edge length
            l = np.sqrt((element.n1.x - element.n2.x)**2 + (element.n1.y - element.n2.y)**2) # length of the edge'

            # integral over N_a*q = l/2 * q
            f_e[0] = l * q * N1(element.n1.x, element.n1.y) # N1 is the basis function for node n1
            f_e[1] = l * q * N2(element.n2.x, element.n2.y) # N2 is the basis function for node n2
            f_e[2] = 0 # N3 is not on the edge            

        elif element.n2 in neumann_nodes and element.n3 in neumann_nodes:
            # edge n2-n3
            l = np.sqrt((element.n2.x - element.n3.x)**2 + (element.n2.y - element.n3.y)**2)
            f_e[0] = 0 # N1 is not on the edge
            f_e[1] = l * q * N2(element.n2.x, element.n2.y) # N2 is the basis function for node n2
            f_e[2] = l * q * N3(element.n3.x, element.n3.y) # N3 is the basis function for node n3
           
        elif element.n3 in neumann_nodes and element.n1 in neumann_nodes:
            # edge n3-n1
            l = np.sqrt((element.n3.x - element.n1.x)**2 + (element.n3.y - element.n1.y)**2)
            f_e[0] = l * q * N1(element.n1.x, element.n1.y) # N1 is the basis function for node n1
            f_e[1] = 0 # N2 is not on the edge
            f_e[2] = l * q * N3(element.n3.x, element.n3.y) # N3 is the basis function for node n3
            
        f[element.n1.id-1] += f_e[0]
        f[element.n2.id-1] += f_e[1]
        f[element.n3.id-1] += f_e[2]

    return -f



def compute_rhs(H, mesh, f, T_dirichlet):
    neumann_nodes = mesh.neumann_nodes
    dirichlet_nodes = mesh.dirichlet_nodes
    neumann_nodes_inside = mesh.neumann_nodes_inside
    free_nodes = np.concatenate([neumann_nodes_inside, neumann_nodes]) 

    
    # get vector P_11...100 for right hand side
    P_1 = f[[node.id -1 for node in free_nodes]]
    
    # get submatrix of H for Neumann nodes (H_1...10,11...100)
    H_neumann = np.delete(H, [node.id-1 for node in dirichlet_nodes], axis=0) # removed rows corresponding to Dirichlet nodes 
    H_neumann = np.delete(H_neumann, [node.id-1 for node in free_nodes], axis=1) # removed columns corresponding to free nodes
    
    T_0 = np.full(len(dirichlet_nodes), T_dirichlet)

    rhs = P_1 - np.matmul(H_neumann, T_0) # right hand side of the system of equations
    
    return rhs

def compute_reaction_forces(H, mesh, T):
    P = np.zeros(len(mesh.nodes)) # initialization of T_1...100
    P = np.matmul(H, T) 

    dirichlet_nodes = mesh.dirichlet_nodes
    return P[[node.id-1 for node in dirichlet_nodes]]# reaction forces at Dirichlet nodes

def temperature_gradient(mesh, T_global):
    """page 124
    values related to the geometrical center of the element

    temp is constant inside the element

    Parameters:
    - mesh object, contains elements and nodes coordinates
    - T: solved global temperature vector

    Returns:
    - grad_T: numpy array [∂T/∂x, ∂T/∂y]

    logic draft:
    loop over elements in the mesh
    obtain the local node ids 
    get local nodal temperatures from t_global
    get b and c coeffs, and area
    build the matrixes and solve the system
    """

    for element in mesh.elements:
        node_ids = element.node_ids()          
        T_local = T_global[node_ids] # shape (3,)

        b = element.b_coeffs() # -> list
        c = element.c_coeffs() # -> list
        A = element.area()

        BC_matrix = np.array([b, c]) # shape: (2, 3)
        grad_T = (1 / (2 * A)) * BC_matrix @ T_local # shape: (2,)

        element.gradient = grad_T


def compute_heat_flux(mesh, k):
    """ q = -k * grad(T) at the element centroids
    k thermal conductivity constant across the element"""
    for element in mesh.elements:
        grad_T = element.gradient
        element.flux = -k * grad_T

def plot_temperature_field(mesh, T):
    """Plot the temperature field T(x, y) as contour plot, with the nodal temperatures as primary
values and (bi-linear) interpolation between the nodes"""

    # to plot triangled mesh, get the x and y nodes and the triangles
    node_x = np.array([node.x for node in mesh.nodes])
    node_y = np.array([node.y for node in mesh.nodes])
    triangles = np.array([
        element.node_ids() for element in mesh.elements
    ])

    triangulation = tri.Triangulation(node_x, node_y, triangles)

    plt.figure(figsize=(8, 6))
    contour = plt.tricontourf(triangulation, T, levels=50, cmap='plasma')
    plt.colorbar(contour, label="Temperature")
    plt.title("Temperature Field $T(x, y)$")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.grid(True)
    plt.show()

def plot_temperature_gradient(mesh):

    centroids_x = np.array([element.centroidX for element in mesh.elements])
    centroids_y = np.array([element.centroidY for element in mesh.elements])
    grad_x = np.array([element.gradient[0]for element in mesh.elements])
    grad_y = np.array([element.gradient[1]for element in mesh.elements])

    # to plot triangled mesh, get the x and y nodes and the triangles
    node_x = np.array([node.x for node in mesh.nodes])
    node_y = np.array([node.y for node in mesh.nodes])
    triangles = np.array([
        element.node_ids() for element in mesh.elements
    ])

    # Plot vector field
    plt.figure(figsize=(8, 6))
    plt.quiver(centroids_x, centroids_y, grad_x, grad_y, angles='xy', scale_units='xy', scale=50, color='blue')
    plt.title("Temperature Gradient Vectors at Element Centroids")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.grid(True)
    plt.show()

    # extra contour for each component of the temp gradient
    triangulation = tri.Triangulation(node_x, node_y, triangles)

    # ∂T/∂x
    plt.figure(figsize=(8, 6))
    tpc = plt.tripcolor(triangulation, facecolors=grad_x, edgecolors='k', shading='flat', cmap='coolwarm')
    plt.colorbar(tpc, label="$\\partial T/\\partial x$")
    plt.title("Constant ∂T/∂x Per Element")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.grid(False)
    plt.show()

    # ∂T/∂y
    plt.figure(figsize=(8, 6))
    tpc = plt.tripcolor(triangulation, facecolors=grad_y, edgecolors='k', shading='flat', cmap='coolwarm')
    plt.colorbar(tpc, label="$\\partial T/\\partial y$")
    plt.title("Constant ∂T/∂y Per Element")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.grid(True)
    plt.show()


def plot_heat_flux(mesh):

    centroids_x = np.array([element.centroidX for element in mesh.elements])
    centroids_y = np.array([element.centroidY for element in mesh.elements])
    flux_x = np.array([element.flux[0]for element in mesh.elements])
    flux_y = np.array([element.flux[1]for element in mesh.elements])

    # to plot triangled mesh, get the x and y nodes and the triangles
    node_x = np.array([node.x for node in mesh.nodes])
    node_y = np.array([node.y for node in mesh.nodes])
    triangles = np.array([
        element.node_ids() for element in mesh.elements
    ])

    # Plot vector field
    plt.figure(figsize=(8, 6))
    plt.quiver(centroids_x, centroids_y, flux_x, flux_y, angles='xy', scale_units='xy', scale=50, color='blue')
    plt.title("Heat Flux Vectors at Element Centroids")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.grid(True)
    plt.show()

    # extra contour for each component of the heat flux
    triangulation = tri.Triangulation(node_x, node_y, triangles)

    # in x
    plt.figure(figsize=(8, 6))
    tpc = plt.tripcolor(triangulation, facecolors=flux_x, edgecolors='k', shading='flat', cmap='coolwarm')
    plt.colorbar(tpc, label="q in x direction")
    plt.title("Heat flux in x Per Element")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.grid(False)
    plt.show()

    # in y
    plt.figure(figsize=(8, 6))
    tpc = plt.tripcolor(triangulation, facecolors=flux_y, edgecolors='k', shading='flat', cmap='coolwarm')
    plt.colorbar(tpc, label="q in y direction")
    plt.title("Heat flux in y Per Element")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.grid(True)
    plt.show()

def plot_temperature_gradients_and_fluxes():
    """Plot the temperature gradients and fluxes at the element centroids as vector plots 
    (and, optionally,their components as contour plots without interpolation, i.e. constant 
    for each element)."""

    pass

def compare_fluxes():
    """Compare the fluxes in elements attached to the boundary y = L to the applied Neumann
BC values (applied via nodal forces P91..100)."""
    pass

