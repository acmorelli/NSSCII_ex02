import numpy as np
from mesh import *
from functions import *
from numpy.linalg import solve
# made by GPT, it doesnt work !!!

def print_HTP(H, T, P):
    print("Stiffness Matrix H:")
    print(np.array_str(H, precision=2, suppress_small=True))
    print("\nTemperature Vector T:")
    print(np.array_str(T, precision=2, suppress_small=True))
    print("\nLoad Vector P:")
    print(np.array_str(P, precision=2, suppress_small=True))

def main_test_v1():
    # Set up the mesh with V1 variation
    mesh = Mesh(N=4, Variation='V1', L=1.0)

    # Override boundary conditions for known analytical solution
    # T = 100 at y=0, T = 0 at y=1, insulated sides (Neumann = 0)
    mesh.T_known = {nid: 100.0 for nid, coord in mesh.nodes.items() if np.isclose(coord[1], 0.0)}
    mesh.T_known.update({nid: 0.0 for nid, coord in mesh.nodes.items() if np.isclose(coord[1], mesh.L)})
    mesh.P_known = {}  # No Neumann fluxes, all insulated

    nodes = mesh.nodes
    elements = mesh.elements
    T_known = mesh.T_known
    P_known = mesh.P_known
    k = mesh.k
    h_z = mesh.hz
    n_nodes = len(nodes)

    # Initialize global matrices
    H = np.zeros((n_nodes, n_nodes))
    P = np.zeros(n_nodes)
    T = np.zeros(n_nodes)

    # Assemble global stiffness matrix
    for elem_id, node_ids in elements.items():
        coords = [nodes[nid] for nid in node_ids]
        He = compute_local_stiffness_matrix(coords, k, h_z)
        for i in range(3):
            for j in range(3):
                H[node_ids[i]-1, node_ids[j]-1] += He[i, j]

    # Apply known Dirichlet BCs
    for nid, val in T_known.items():
        T[nid - 1] = val

    known_ids = list(T_known.keys())
    free_ids = [i for i in range(1, n_nodes+1) if i not in known_ids]

    # Solve linear system for unknown temperatures
    H_ff = H[np.ix_([i-1 for i in free_ids], [i-1 for i in free_ids])]
    H_fk = H[np.ix_([i-1 for i in free_ids], [i-1 for i in known_ids])]
    P_f = P[[i-1 for i in free_ids]] - H_fk @ T[[i-1 for i in known_ids]]

    T_f = solve(H_ff, P_f)
    for idx, nid in enumerate(free_ids):
        T[nid - 1] = T_f[idx]

    # Compute reaction forces at known BC nodes
    P = H @ T

    # Output results
    print_HTP(H, T, P)

    total_reaction = sum(P[i-1] for i in known_ids)
    print(f"\nTotal Reaction Heat: {total_reaction:.2f} (Expected: ~100.0)")

if __name__ == "__main__":
    main_test_v1()
