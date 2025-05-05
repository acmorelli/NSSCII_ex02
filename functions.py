from min_FEM import *
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri
import os

def compute_local_stiffness_matrix(k,element_id, mesh, hz,Variation=None, id_c=None, ce=None):

    # k: heat conductivity coefficient
    # tri_id: id of the triangle in the mesh
    # mesh: mesh object

    H_e = np.zeros((3, 3)) # local stiffness matrix


    # get element
    element = mesh.elements[element_id-1]
    # coefficients
    b = np.array(element.b_coeffs())
    c = np.array(element.c_coeffs())
    # Area of the element
    A = element.area()
    if Variation == 'V4a':
        if element.id in id_c:
            element.k=ce*element.k
        else:
            pass
    elif Variation == 'V4b':
        if element.id in id_c:
            element.k=element.k/ce
        else:
            pass
    else:
        pass
    # see Zienkewicz, page 120 and 125
    for i in range(3):
        for j in range(3):
            H_e[i][j] = (element.k* hz / (4*A)) * (b[i]*b[j] + c[i]*c[j]) 

    return H_e

def compute_free_nodes(H, mesh):
    dirichlet_nodes = mesh.dirichlet_nodes
    # matrix without Dirichlet nodes (free nodes)
    H_free = np.delete(H, [node.id-1 for node in dirichlet_nodes], axis=0)
    H_free = np.delete(H_free, [node.id-1 for node in dirichlet_nodes], axis=1)

    return H_free

def compute_load_vector(H, mesh, q, hz):
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
            print('edge n1-n2', element.id, element.n1, element.n2)
            # edge length
            l = np.sqrt((element.n1.x - element.n2.x)**2 + (element.n1.y - element.n2.y)**2) # length of the edge'

            # integral over N_a*q = l/2 * q
            # f_e[0] = l * q * N1(element.n1.x, element.n1.y) # N1 is the basis function for node n1
            # f_e[1] = l * q * N2(element.n2.x, element.n2.y) # N2 is the basis function for node n2
            f_e[0] =0.5 * l * q * hz #* N1(element.n1.x, element.n1.y)
            f_e[1] = 0.5*l * q * hz #* N2(element.n2.x, element.n2.y)
            f_e[2] = 0 # N3 is not on the edge

        elif element.n2 in neumann_nodes and element.n3 in neumann_nodes:
            # edge n2-n3
            print('edge n2-n3', element.id, element.n2, element.n3)
            l = np.sqrt((element.n2.x - element.n3.x)**2 + (element.n2.y - element.n3.y)**2)
            f_e[0] = 0 # N1 is not on the edge
            f_e[1] = 0.5 * l * q * hz #N2(element.n2.x, element.n2.y) * hz# N2 is the basis function for node n2
            f_e[2] = 0.5 * l * q * hz #N3(element.n3.x, element.n3.y) * hz# N3 is the basis function for node n3

        elif element.n3 in neumann_nodes and element.n1 in neumann_nodes:
            print('edge n3-n1', element.id, element.n3, element.n1)

            # edge n3-n1
            l = np.sqrt((element.n3.x - element.n1.x)**2 + (element.n3.y - element.n1.y)**2)
            f_e[0] = 0.5 * l * q * hz #N1(element.n1.x, element.n1.y)* hz # N1 is the basis function for node n1
            f_e[1] = 0 # N2 is not on the edge
            f_e[2] = 0.5 * l * q * hz #N3(element.n3.x, element.n3.y)* hz # N3 is the basis function for node n3

        f[element.n1.id-1] += f_e[0]
        f[element.n2.id-1] += f_e[1]
        f[element.n3.id-1] += f_e[2]
    return -f



def compute_rhs(H, mesh, f, T_dirichlet):
    neumann_nodes = mesh.neumann_nodes
    dirichlet_nodes = mesh.dirichlet_nodes
    neumann_nodes_inside = mesh.neumann_nodes_inside
    # free_nodes = np.concatenate([neumann_nodes_inside, neumann_nodes])
    free_nodes = np.array([node for node in mesh.nodes if node not in dirichlet_nodes])
    
    # get vector P_11...100 for right hand side
    P_1 = f[[node.id -1 for node in free_nodes]]
    print("f", f)
    print('P_1', P_1)
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

def compute_temperature_gradient(mesh, T_global):
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
        grad_T = (1 / (2*A)) * BC_matrix @ T_local # shape: (2,)
        # not divide by two because the coefficients are derived assuming area-normalized shape functions
        # print('grad_T', grad_T)
        element.gradient = grad_T


def compute_heat_flux(mesh):
    """ q = -k * grad(T) at the element centroids
    k thermal conductivity constant across the element"""
    for element in mesh.elements:
        #grad_T = element.gradient
        element.flux = -element.k * element.gradient

def plot_temperature_field(mesh, T, variation):
    """Plot the temperature field T(x, y) as contour plot, with the nodal temperatures as primary
values and (bi-linear) interpolation between the nodes"""

    # to plot triangled mesh, get the x and y nodes and the triangles
    node_x = np.array([node.x for node in mesh.nodes]) # global
    node_y = np.array([node.y for node in mesh.nodes])
    triangles = np.array([
        element.node_ids() for element in mesh.elements
    ])

    triangulation = tri.Triangulation(node_x, node_y, triangles)

    fig=plt.figure(figsize=(8, 6))
    contour = plt.tricontourf(triangulation, T, levels=50, cmap='plasma')
    plt.colorbar(contour, label="Temperature")
    plt.title("Temperature Field $T(x, y)$")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.grid(False)
    plt.show()

    folder = f"./{variation}_plots"
    os.makedirs(folder, exist_ok=True)
    fig.savefig(os.path.join(folder, "temp_field_plot.png"), bbox_inches='tight')
    plt.close(fig)

def plot_temperature_gradient(mesh, variation):
    folder = f"./{variation}_plots"
    os.makedirs(folder, exist_ok=True)

    # get centroids, gradients and triangle
    centroids_x = np.array([element.centroidX for element in mesh.elements])
    centroids_y = np.array([element.centroidY for element in mesh.elements])
    grad_x = np.array([element.gradient[0]for element in mesh.elements])
    grad_y = np.array([element.gradient[1]for element in mesh.elements])
    print('grad y', grad_y)
    node_x = np.array([node.x for node in mesh.nodes])
    node_y = np.array([node.y for node in mesh.nodes])
    triangles = np.array([
        element.node_ids() for element in mesh.elements
    ])

    # plots
    fig=plt.figure(figsize=(8, 6))
    plt.quiver(centroids_x, centroids_y, grad_x, grad_y, angles='xy', scale_units='xy', scale=100000, color='blue')
    plt.title("Temperature Gradient Vectors at Elements' Centroids")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.grid(False)
    plt.show()    
    
    fig.savefig(os.path.join(folder, "temp_grad_vector_plot.png"), bbox_inches='tight')
    plt.close(fig)

    # extra contour for each component of the temp gradient
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    triangulation = tri.Triangulation(node_x, node_y, triangles)
   # ∂T/∂x
    tpc1 = axs[0].tripcolor(triangulation, facecolors=np.round(grad_x, 7), edgecolors='k',
                         shading='flat', cmap='viridis')
    fig.colorbar(tpc1, ax=axs[0], label="$\\partial T/\\partial x$")
    axs[0].set_title("∂T/∂x Per Element")
    axs[0].set_xlabel("X")
    axs[0].set_ylabel("Y")
    axs[0].axis("equal")

    # ∂T/∂y
    tpc2 = axs[1].tripcolor(triangulation, facecolors=np.round(grad_y, 7), edgecolors='k',
                         shading='flat', cmap='plasma')
    fig.colorbar(tpc2, ax=axs[1], label="$\\partial T/\\partial y$")
    axs[1].set_title("∂T/∂y Per Element")
    axs[1].set_xlabel("X")
    axs[1].set_ylabel("Y")
    axs[1].axis("equal")
    plt.show()

    fig.savefig(os.path.join(folder, "components_temp_grad.png"), bbox_inches='tight')


def plot_heat_flux(mesh, variation):
    folder = f"./{variation}_plots"
    os.makedirs(folder, exist_ok=True)

    centroids_x = np.array([element.centroidX for element in mesh.elements])
    centroids_y = np.array([element.centroidY for element in mesh.elements])
    flux_x = np.array([element.flux[0]for element in mesh.elements])
    flux_y = np.array([element.flux[1]for element in mesh.elements])
    print('flux y', flux_y)
    # to plot triangled mesh, get the x and y nodes and the triangles
    node_x = np.array([node.x for node in mesh.nodes])
    node_y = np.array([node.y for node in mesh.nodes])
    triangles = np.array([
        element.node_ids() for element in mesh.elements
    ])

    # Plot vector field
    fig=plt.figure(figsize=(8, 6))
    plt.quiver(centroids_x, centroids_y, flux_x, flux_y, angles='xy', scale_units='xy', scale=100000000, color='blue')
    plt.title("Heat Flux Vectors at Element Centroids")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.grid(False)
    plt.show()

    fig.savefig(os.path.join(folder, "flux_vector_plot.png"), bbox_inches='tight')
    plt.close(fig)

    # extra contour for each component of the heat flux
    triangulation = tri.Triangulation(node_x, node_y, triangles)

    fig, axs = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

    # in x
    rounded_flux_x = np.array(np.round(flux_x, decimals=7))
    tpc_x = axs[0].tripcolor(triangulation, facecolors=rounded_flux_x, edgecolors='k',
                             shading='flat', cmap='viridis')
    fig.colorbar(tpc_x, ax=axs[0], label="q in x direction")
    axs[0].set_title("Heat Flux in X Per Element")
    axs[0].set_xlabel("X")
    axs[0].set_ylabel("Y")
    axs[0].axis('equal')
    axs[0].grid(False)

    # in y
    rounded_flux_y= np.array(np.round(flux_y, decimals=7))
    tpc_y = axs[1].tripcolor(triangulation, facecolors=rounded_flux_y, edgecolors='k',
                             shading='flat', cmap='plasma')
    fig.colorbar(tpc_y, ax=axs[1], label="q in y direction")
    axs[1].set_title("Heat Flux in Y Per Element")
    axs[1].set_xlabel("X")
    axs[1].set_ylabel("Y")
    axs[1].axis('equal')
    axs[1].grid(False)

    plt.show()
    fig.savefig(os.path.join(folder, "flux_components_plot.png"), bbox_inches='tight')
    plt.close(fig)

def plot_temperature_gradients_and_fluxes():
    """Plot the temperature gradients and fluxes at the element centroids as vector plots 
    (and, optionally,their components as contour plots without interpolation, i.e. constant 
    for each element)."""

    pass

def compare_fluxes():
    """Compare the fluxes in elements attached to the boundary y = L to the applied Neumann
BC values (applied via nodal forces P91..100)."""
    pass

