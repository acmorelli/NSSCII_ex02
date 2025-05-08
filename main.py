"""observations version 1 - squared mesh:
• The temperature gradient must be constant in the entire model and equal to the overall
gradient ΔT/Δy.
• The flux must be constant in the entire model and must be equal to the applied flux. The
latter must be equal to the sum of the applied nodal “forces” divided by the corresponding
area.
• The temperature gradient and the flux, together with eqn. (2) must yield the input for
the conductivity k.

SI-units to be used:
    + T in K
    + L in m
    + k in W/(mK)
    + q in W/m^2 - ad Neumann
    + P in W - ad nodal forces

# Length in x- and y-direction.
L = 0.1

# Thickness (z-direction).
hz = 0.01

# Thermal conductivity (k=k_xx=k_yy, k_xy = 0.).
k = 373.

# Factor c for modifying thermal conductivity k for
# elements in elements_to_be_modified.
c = 10.

# Elements to be modified.
elements_to_be_modified = [
                          62-70,
                          83-86,
                          98-103,
                          111-122
                          ]

# Boundary conditions.
q(y=0) = 1000000.
T(y=L) = 313.
"""

# import local files here
from functions import *
from mesh import Mesh
import numpy as np
import matplotlib.pyplot as plt
from print_HTP_2025 import print_HTP


def main():

    k = 373  #  W/mKn
    # convert to W/K
    # k = k * 1e3 # W/K

    L = 0.1  # m (length of squared domain (V0))
    N = 100  # Number of nodes in x=0 (V0)
    hz = 0.01  # m (thickness in z-direction)
    q_neumann = 1000000 * hz  # m* w/ m2 (Flux across the Neumann boundary)
    y_neumann = 0.0  # y coordinate of the Neumann boundary
    T_dirichlet = 313.0  # K (Dirichlet bc)
    y_dirichlet = L  # y coordinate of the Dirichlet boundary
    Variation = "V4b"  # 'V0', 'V1', 'V2', 'V3', 'V4a', 'V4b' (for the mesh)

    """
    # # # debug
    L=1
    k=1
    q_neumann=1
    y_neumann= L
    T_dirichlet = 1
    y_dirichlet = 0.0
    """

    # # # # variation 4 # # # #
    id1 = np.arange(62, 71)
    id2 = np.arange(83, 87)
    id3 = np.arange(98, 104)
    id4 = np.arange(111, 123)
    id_c = np.concatenate((id1, id2, id3, id4))
    ce = 10.0  # factor for modifying thermal conductivity k
    # # # # # # # # # # # # # # 

    mesh = Mesh(N, Variation, L, k, y_neumann, y_dirichlet, hz)
    H = np.zeros((N, N))  # Global stiffness matrix initialization
    for element_id in range(1, len(mesh.elements) + 1):

        # Compute local stiffness matrix for each element
        H_e = compute_local_stiffness_matrix(
            k, element_id, mesh, hz, Variation, id_c, ce
        )

        # Assemble global stiffness matrix
        element_nodes_ids = [
            mesh.elements[element_id - 1].n1.id,
            mesh.elements[element_id - 1].n2.id,
            mesh.elements[element_id - 1].n3.id,
        ]
        for row_e in range(H_e.shape[0]):
            global_row_idx = element_nodes_ids[row_e] - 1
            for col_e in range(H_e.shape[1]):
                global_col_idx = element_nodes_ids[col_e] - 1
                H[global_row_idx][global_col_idx] += H_e[row_e][col_e]

    H_free = compute_free_nodes(H, mesh)

    f = compute_load_vector(H, mesh, q_neumann, hz)

    rhs = compute_rhs(H, mesh, f, T_dirichlet)

    # solve the system of equations
    T_free = np.linalg.solve(H_free, rhs)

    # set up complete solution vector
    T = np.zeros(len(mesh.nodes))

    free_nodes = np.array(
        [node for node in mesh.nodes if node not in mesh.dirichlet_nodes]
    )
    for i, node in enumerate(free_nodes):
        T[node.id - 1] = T_free[i]
        #print("node", node, "t_free[i]", T_free[i])
    for node in mesh.dirichlet_nodes:
        T[node.id - 1] = T_dirichlet
    #print("T:", T)

    reaction_forces = compute_reaction_forces(H, mesh, T)

    # add reaction forces to force vector
    counter = 0
    for node in mesh.nodes:
        if node in mesh.dirichlet_nodes:
            f[node.id - 1] += reaction_forces[counter]
            counter += 1

    # post processing
    compute_temperature_gradient(mesh, T)
    compute_heat_flux(mesh)

    # plotting
    plot_temperature_field(mesh, T, Variation)
    plot_temperature_gradient(mesh, Variation)
    plot_heat_flux(mesh, Variation)

    # save results to file
    filename = "output_" + Variation + ".txt"
    print_HTP(H, np.reshape(T, (-1, 1)), np.reshape(f, (-1, 1)), filename)
    print(f"Results saved to {filename}")


if __name__ == "__main__":
    main()
