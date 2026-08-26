import yaml
from pathlib import Path


def _get_project_root() -> Path:
    return Path(__file__).resolve().parents[3]


def load_annotation_config(config_path=None) -> dict:
    if config_path is None:
        project_root = _get_project_root()
        config_path = project_root / "config" / "annotation.yaml"

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    _validate_config(cfg)
    return cfg


def _validate_config(cfg: dict):
    required_sections = ["paths", "filter", "annotation", "logging"]
    for section in required_sections:
        if section not in cfg:
            raise KeyError(f"Missing required config section: '{section}'")

    paths = cfg["paths"]
    for key in ["csv_path", "cif_dir", "output_db_dir"]:
        if key not in paths:
            raise KeyError(f"Missing required path config: 'paths.{key}'")

    filt = cfg["filter"]
    for key in ["holo_only", "allowed_agtypes", "drop_ion_only",
                 "min_ag_seq_length", "allowed_antibody_types"]:
        if key not in filt:
            raise KeyError(f"Missing required filter config: 'filter.{key}'")

    ann = cfg["annotation"]
    for key in ["contact_cutoff", "rsa_threshold", "sasa_point_number",
                 "min_chain_residues"]:
        if key not in ann:
            raise KeyError(f"Missing required annotation config: 'annotation.{key}'")
