import yaml
from pathlib import Path

import pandas as pd


def load_filter_config(config_path=None):
    if config_path is None:
        project_root = Path(__file__).resolve().parents[2]
        config_path = project_root / "config" / "sabdab_filter.yaml"

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg


def _parse_types(types_str):
    if pd.isna(types_str):
        return []
    return [t for t in types_str.split("/") if t]


def _has_only_allowed_types(types, allowed_set):
    return len(types) > 0 and all(t in allowed_set for t in types)


def _is_only_ions(types):
    return len(types) > 0 and all(t == "ION" for t in types)


def filter_holo(df):
    df = df.loc[df["holo"] == True].copy()
    return df


def filter_agtypes(df, config):
    allowed_set = set(config["always_keep_agtypes"] + config.get("additional_agtypes", []))

    mask_allowed = df["agtypes"].apply(lambda s: _has_only_allowed_types(_parse_types(s), allowed_set))
    mask_not_ions = df["agtypes"].apply(lambda s: not _is_only_ions(_parse_types(s)))

    df = df.loc[mask_allowed & mask_not_ions].copy()
    return df


def filter_min_sequence_length(df, config):
    min_len = config["min_sequence_length"]

    mask = df["agresolvedseqs"].apply(
        lambda s: any(len(chain) > min_len for chain in s.split("/"))
        if not pd.isna(s) else False
    )
    df = df.loc[mask].copy()
    return df


def filter_antibody_types(df, config):
    allowed = set(config["allowed_antibody_types"])
    df = df.loc[df["type"].isin(allowed)].copy()
    return df


def drop_columns(df, config, phase):
    cols = config.get(f"columns_to_drop_{phase}", [])
    df = df.drop(columns=[c for c in cols if c in df.columns], axis=1)
    return df


def apply_sabdab_filter(df, config=None, config_path=None):
    if config is None:
        config = load_filter_config(config_path)

    df = drop_columns(df, config, "initial")
    df = filter_holo(df)
    df = filter_agtypes(df, config)
    df = filter_min_sequence_length(df, config)
    df = filter_antibody_types(df, config)
    df = drop_columns(df, config, "after")

    return df
