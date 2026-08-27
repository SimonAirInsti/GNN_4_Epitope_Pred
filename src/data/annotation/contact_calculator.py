import numpy as np
from scipy.spatial import cKDTree


def _build_tree(coords_dict):
    """Flatten {chain: {resseq: [coords]}} into (array, atom→(chain, resseq) mapping)."""
    atoms = []
    mapping = []
    for chain, res_dict in coords_dict.items():
        for resseq, atom_list in res_dict.items():
            for atom in atom_list:
                atoms.append(atom)
                mapping.append((chain, resseq))
    return np.array(atoms, dtype=np.float64), mapping


def compute_distances(query_coords, target_coords, cutoff=4.5):
    """
    Compute all pairwise distances between query and target residues.

    Builds a single KDTree from target atoms and queries every query residue
    against it. Returns distances, contact pairs, and per-side residue distances
    in one pass.

    Args:
        query_coords: {chain: {resseq: [coords]}} residues to measure from.
        target_coords: {chain: {resseq: [coords]}} residues to measure against.
        cutoff: Distance threshold (Angstroms) for contact pairs.

    Returns:
        dict with:
            - query_dists: {chain: {resseq: min_distance}} for query residues.
            - contact_pairs: list of (q_chain, q_resseq, t_chain, t_resseq, distance).
    """
    target_array, target_mapping = _build_tree(target_coords)
    tree = cKDTree(target_array)

    query_dists = {}
    contact_pairs = []

    for q_chain, res_dict in query_coords.items():
        query_dists[q_chain] = {}
        for q_resseq, atom_list in res_dict.items():
            q_array = np.array(atom_list, dtype=np.float64)
            dists, indices = tree.query(q_array, k=1)
            min_dist = float(np.min(dists))
            query_dists[q_chain][q_resseq] = min_dist

            if min_dist <= cutoff:
                best_idx = int(np.argmin(dists))
                t_chain, t_resseq = target_mapping[indices[best_idx]]
                contact_pairs.append(
                    (q_chain, q_resseq, t_chain, t_resseq, min_dist)
                )

    return {
        "query_dists":    query_dists,
        "contact_pairs":  contact_pairs,
    }


def compute_contacts(ag_coords, ab_coords):
    """Compute min distance from each antigen residue to any antibody atom."""
    return compute_distances(ag_coords, ab_coords)["query_dists"]


def compute_paratope(ab_coords, ag_coords):
    """Compute min distance from each antibody residue to any antigen atom."""
    return compute_distances(ab_coords, ag_coords)["query_dists"]
