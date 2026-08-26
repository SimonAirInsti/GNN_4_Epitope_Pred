import os
from collections import defaultdict
from Bio import PDB


def load_complex_coords(row, cif_dir, min_chain_residues=5, logger=None):
    """
    Parse a CIF file and extract heavy-atom coordinates + sequences for
    antibody and antigen chains.

    Drops antigen chains with fewer than min_chain_residues residues
    (typically noise ions).

    Args:
        row: DataFrame row with INSTANCE, Hchain, Lchain, agchains columns.
        cif_dir: Directory containing .cif files.
        min_chain_residues: Minimum residues for an antigen chain to be kept.
        logger: Optional logging.Logger instance.

    Returns:
        dict with cif_path, instance, antibody/antigen_chains,
                antibody/antigen_coords, antibody/antigen_sequences.
    """
    instance = row["INSTANCE"]
    cif_path = os.path.join(cif_dir, f"{instance}.cif")

    if not os.path.isfile(cif_path):
        raise FileNotFoundError(f"CIF not found: {cif_path}")

    # --- Parse CIF -> heavy atoms per residue ---
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("complex", cif_path)
    model = structure[0]

    all_coords = defaultdict(lambda: defaultdict(list))

    for residue in model.get_residues():
        chain_id = residue.get_parent().get_id()
        resseq = residue.get_id()[1]
        for atom in residue:
            if atom.element.upper() == "H":
                continue
            all_coords[chain_id][resseq].append(tuple(atom.get_coord()[:3]))

    all_coords = dict(all_coords)

    # --- Extract amino acid sequences ---
    ppb = PDB.PPBuilder()
    all_sequences = {}
    for mdl in structure:
        for chain in mdl:
            polypeptides = ppb.build_peptides(chain)
            if polypeptides:
                seq = "".join(str(pp.get_sequence()) for pp in polypeptides)
                all_sequences[chain.id] = seq
            else:
                all_sequences[chain.id] = ""

    # --- Identify antibody chains ---
    h_chain = str(row["Hchain"]).strip()
    l_chain = str(row["Lchain"]).strip()

    antibody_chains = set()
    if h_chain and h_chain not in ("nan", "None", ""):
        antibody_chains.add(h_chain)
    if l_chain and l_chain not in ("nan", "None", "") and l_chain != h_chain:
        antibody_chains.add(l_chain)

    # --- Identify antigen chains ---
    ag_chain_str = str(row["agchains"]).strip()
    antigen_chains = set(c.strip() for c in ag_chain_str.split("/") if c.strip())

    # --- Filter to chains that exist in CIF ---
    cif_chains = set(all_coords.keys())
    antibody_chains &= cif_chains
    antigen_chains &= cif_chains

    if antibody_chains & antigen_chains:
        overlap = antibody_chains & antigen_chains
        if logger:
            logger.warning(f"Chain overlap Ab/Ag in {instance}: {overlap}")

    if not antibody_chains:
        raise ValueError(
            f"{instance}: no antibody chains found "
            f"({h_chain}/{l_chain}) in CIF"
        )
    if not antigen_chains:
        raise ValueError(
            f"{instance}: no antigen chains found "
            f"({ag_chain_str}) in CIF"
        )

    # --- Filter antigen chains by minimum residue count ---
    dropped_chains = set()
    for ch in list(antigen_chains):
        n_res = len(all_coords.get(ch, {}))
        if n_res < min_chain_residues:
            dropped_chains.add(ch)
            antigen_chains.discard(ch)
            if logger:
                logger.warning(
                    f"{instance}: antigen chain {ch} dropped "
                    f"({n_res} residues < min_chain_residues={min_chain_residues})"
                )

    if dropped_chains:
        dropped_msg = ", ".join(sorted(dropped_chains))
        if logger:
            logger.info(
                f"{instance}: dropped antigen chains: {dropped_msg}"
            )

    if not antigen_chains:
        raise ValueError(
            f"{instance}: all antigen chains dropped below "
            f"min_chain_residues threshold ({ag_chain_str})"
        )

    return {
        "cif_path":           cif_path,
        "instance":           instance,
        "antibody_chains":    antibody_chains,
        "antigen_chains":     antigen_chains,
        "antibody_coords":    {ch: all_coords[ch] for ch in antibody_chains},
        "antigen_coords":     {ch: all_coords[ch] for ch in antigen_chains},
        "antibody_sequences": {ch: all_sequences.get(ch, "") for ch in antibody_chains},
        "antigen_sequences":  {ch: all_sequences.get(ch, "") for ch in antigen_chains},
        "dropped_chains":     dropped_chains,
    }
