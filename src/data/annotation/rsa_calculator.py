import numpy as np
from collections import defaultdict
from biotite.structure.io import load_structure
from biotite.structure import sasa


MAX_SSA = {
    "ALA": 107.95, "ARG": 248.34, "ASN": 157.99, "ASP": 153.80,
    "CYS": 140.38, "GLN": 198.89, "GLU": 194.70, "GLY": 84.21,
    "HIS": 189.34, "ILE": 167.64, "LEU": 169.81, "LYS": 203.51,
    "MET": 178.27, "PHE": 197.11, "PRO": 137.28, "SER": 124.23,
    "THR": 140.05, "TRP": 242.62, "TYR": 223.70, "VAL": 144.86,
}

def compute_rsa(cif_path, point_number=500):
    """
    Compute RSA for each residue using biotite's Shrake-Rupley algorithm.

    Args:
        cif_path: Path to the CIF file.
        point_number: Number of sampling points for SASA calculation.

    Returns:
        dict {chain_id: {resseq: rsa_value}}
    """
    array = load_structure(cif_path)
    atom_sasa = sasa(array, point_number=point_number)

    residue_sasa = defaultdict(lambda: [0.0, None])

    for j in range(len(array)):
        sasa_j = atom_sasa[j]
        if not np.isfinite(sasa_j):
            continue
        chain_id = array.chain_id[j]
        resseq = array.res_id[j]
        resname = array.res_name[j]
        key = (chain_id, resseq)
        residue_sasa[key][0] += sasa_j
        residue_sasa[key][1] = resname

    rsa = defaultdict(dict)
    for (chain_id, resseq), (sasa_real, resname) in residue_sasa.items():
        max_area = MAX_SSA.get(resname, 150.0)
        if max_area > 0:
            rsa[chain_id][resseq] = sasa_real / max_area
        else:
            rsa[chain_id][resseq] = 0.0

    return dict(rsa)
