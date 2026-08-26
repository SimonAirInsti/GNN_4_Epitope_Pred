import pandas as pd


def apply_filters(df, filter_cfg):
    """
    Filter the SAbDab2 DataFrame to keep only valid antibody-antigen complexes.

    Applied before CIF parsing to avoid loading structures that will be discarded.

    Filter chain (applied in order):
        1. Drop depositor columns (PDBdepo, SABDABdepo, SABDABupdate)
        2. Keep only holo conformations (holo == True)
        3. Keep only complexes with allowed antigen types
        4. Drop complexes where agtypes is exclusively ION
        5. Keep complexes with at least one antigen chain > min_ag_seq_length
        6. Keep complexes with allowed antibody types (FAB, FV)
        7. Drop post-filter columns (resolution, holo, method, etc.)

    Args:
        df: Raw SAbDab2 DataFrame.
        filter_cfg: The 'filter' section from annotation.yaml.

    Returns:
        Filtered DataFrame.
    """
    # Drop depositor columns no longer needed
    df = df.drop(
        columns=["PDBdepo", "SABDABdepo", "SABDABupdate"],
        errors="ignore",
    )

    # Holo conformations only (non-apo)
    if filter_cfg.get("holo_only", True):
        df = df.loc[df["holo"] == True].copy()

    # Allowed antigen types: each "/"-separated type must be in the allowed set
    allowed = set(filter_cfg.get("allowed_agtypes", ["PROTEIN", "PEPTIDE", "ION"]))
    mask = df["agtypes"].apply(
        lambda s: all(t in allowed for t in str(s).split("/"))
        if pd.notna(s) else False
    )
    df = df.loc[mask].copy()

    # Drop complexes where agtypes is exclusively ION
    if filter_cfg.get("drop_ion_only", True):
        ion_mask = df["agtypes"].apply(
            lambda s: all(
                t.strip() == "ION" for t in str(s).split("/") if t.strip()
            ) if pd.notna(s) else False
        )
        df = df.loc[~ion_mask].copy()

    # At least one antigen chain must exceed min_ag_seq_length
    min_len = filter_cfg.get("min_ag_seq_length", 25)
    seq_mask = df["agresolvedseqs"].apply(
        lambda s: any(len(c) > min_len for c in str(s).split("/"))
        if pd.notna(s) else False
    )
    df = df.loc[seq_mask].copy()

    # Allowed antibody types (FAB, FV are classical paired-chain formats)
    allowed_types = set(filter_cfg.get("allowed_antibody_types", ["FAB", "FV"]))
    df = df.loc[df["type"].isin(allowed_types)].copy()

    # Drop columns no longer needed after filtering
    df = df.drop(
        columns=["resolution", "holo", "method", "species", "construct", "type"],
        errors="ignore",
    )

    return df


def load_csv(csv_path, filter_cfg=None):
    """
    Load SAbDab2 CSV and apply filters.

    Convenience function that combines pd.read_csv with apply_filters.

    Args:
        csv_path: Path to the SAbDab2 split CSV.
        filter_cfg: The 'filter' section from annotation.yaml.

    Returns:
        Filtered DataFrame.
    """
    df = pd.read_csv(csv_path)
    if filter_cfg is not None:
        df = apply_filters(df, filter_cfg)
    return df
