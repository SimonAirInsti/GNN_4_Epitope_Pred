import json

def label_epitope(ag_coords, ag_dists, rsa, cutoff=4.5, rsa_threshold=0.0):
    """
    Label residues as epitope: contact (<= cutoff) AND RSA > threshold.

    Args:
        ag_coords: {chain: {resseq: [coords]}} antigen coords.
        ag_dists: {chain: {resseq: min_distance}} pre-computed Ag→Ab distances.
        rsa: {chain: {resseq: rsa_value}} RSA values.
        cutoff: Contact distance threshold (Angstroms).
        rsa_threshold: Minimum RSA for epitope classification.

    Returns:
        {chain: {resseq: {contact, epitope, distance, rsa}}}
    """
    labels = {}
    for chain in ag_coords:
        labels[chain] = {}
        chain_rsa = rsa.get(chain, {})
        for resseq in ag_coords[chain]:
            d = ag_dists.get(chain, {}).get(resseq, float("inf"))
            r = chain_rsa.get(resseq, 0.0)
            is_contact = d <= cutoff
            is_epitope = is_contact and r > rsa_threshold
            labels[chain][resseq] = {
                "contact": is_contact,
                "epitope": is_epitope,
                "distance": d,
                "rsa": r,
            }
    return labels


def format_epitope_labels(epitope_labels):
    """
    Format epitope labels as 'resseq-chain, resseq-chain, ...'.

    Args:
        epitope_labels: {chain: {resseq: {epitope: bool, ...}}}

    Returns:
        Comma-separated string of epitope residues.
    """
    parts = []
    for chain in sorted(epitope_labels.keys()):
        for resseq in sorted(epitope_labels[chain].keys()):
            if epitope_labels[chain][resseq]["epitope"]:
                parts.append(f"{resseq}-{chain}")
    return ", ".join(parts)


def count_epitope_labels(epitope_labels):
    """Count total epitope residues across all chains."""
    return sum(
        1 for chain in epitope_labels.values()
        for res in chain.values()
        if res["epitope"]
    )


def format_paratope_labels(paratope_dists, cutoff=4.5):
    """
    Format paratope labels as 'resseq-chain, resseq-chain, ...'.

    Args:
        paratope_dists: {chain: {resseq: min_distance}} paratope distances.
        cutoff: Contact distance threshold (Angstroms).

    Returns:
        Comma-separated string of paratope residues.
    """
    parts = []
    for chain in sorted(paratope_dists.keys()):
        for resseq in sorted(paratope_dists[chain].keys()):
            if paratope_dists[chain][resseq] <= cutoff:
                parts.append(f"{resseq}-{chain}")
    return ", ".join(parts)


def count_paratope_labels(paratope_dists, cutoff=4.5):
    """Count total paratope residues across all chains."""
    return sum(
        1 for chain in paratope_dists.values()
        for d in chain.values()
        if d <= cutoff
    )


def serialize_sequences(seq_dict):
    """Convert {chain: sequence} dict to JSON string for DB storage."""
    return json.dumps(seq_dict)
