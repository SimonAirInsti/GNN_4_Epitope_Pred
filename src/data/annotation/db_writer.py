import os
import sqlite3
import json
from datetime import datetime


EPITOME_COLUMNS = [
    "INSTANCE", "PDB_ID", "SABDAB_ID", "HEAVY_ID", "LIGHT_ID",
    "Hchain", "Lchain", "agchains", "agtypes", "agresolvedseqs",
    "agexpectedseqs",
    "cdrh3_cluster", "cdrh123_cluster", "cdrl123_cluster",
    "ab_cluster", "agclusters", "ab_ag_cluster",
    "ab_ag_split",
    "antigen_sequences", "antibody_sequences",
    "epitope_labels", "paratope_labels",
    "n_epitope_residues", "n_paratope_residues",
    "n_contact_pairs", "n_ag_residues",
    "contact_cutoff", "rsa_threshold",
    "status", "error_message",
]


def _table_name_from_config(ann_cfg):
    return f"epitopes_c{ann_cfg['contact_cutoff']}_rsa{ann_cfg['rsa_threshold']}"


def create_database(db_path, ann_cfg):
    """
    Create SQLite database with epitope table and run_metadata table.

    Args:
        db_path: Path to output .db file.
        ann_cfg: The 'annotation' section of the config.

    Returns:
        table_name: Name of the created epitope table.
    """
    table_name = _table_name_from_config(ann_cfg)
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    col_defs = []
    for col in EPITOME_COLUMNS:
        col_type = _infer_type(col)
        if col == "INSTANCE":
            col_defs.append(f"{col} TEXT PRIMARY KEY")
        else:
            col_defs.append(f"{col} {col_type}")

    cur.execute(
        f"CREATE TABLE IF NOT EXISTS [{table_name}] ({', '.join(col_defs)})"
    )

    cur.execute(
        "CREATE TABLE IF NOT EXISTS run_metadata ("
        "run_timestamp TEXT, "
        "config_json TEXT, "
        "table_name TEXT, "
        "total_processed INTEGER, "
        "n_success INTEGER, "
        "n_error INTEGER, "
        "n_skipped INTEGER)"
    )

    conn.commit()
    conn.close()
    return table_name


def _infer_type(col):
    text_cols = [
        "INSTANCE", "PDB_ID", "SABDAB_ID", "HEAVY_ID", "LIGHT_ID",
        "Hchain", "Lchain", "agchains", "agtypes", "agresolvedseqs",
        "agexpectedseqs",
        "cdrh3_cluster", "cdrh123_cluster", "cdrl123_cluster",
        "ab_cluster", "agclusters", "ab_ag_cluster",
        "ab_ag_split",
        "antigen_sequences", "antibody_sequences",
        "epitope_labels", "paratope_labels",
        "status", "error_message",
    ]
    int_cols = [
        "n_epitope_residues", "n_paratope_residues",
        "n_contact_pairs", "n_ag_residues",
    ]
    if col in text_cols:
        return "TEXT"
    if col in int_cols:
        return "INTEGER"
    return "REAL"


def insert_rows(db_path, table_name, rows, full_config):
    """
    Batch-insert annotated rows into the epitope table and write metadata.

    Args:
        db_path: Path to .db file.
        table_name: Epitope table name.
        rows: List of dicts with EPITOME_COLUMNS keys.
        full_config: Full config dict (for metadata snapshot).
    """
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    placeholders = ", ".join(["?" for _ in EPITOME_COLUMNS])
    col_names = ", ".join([f"[{c}]" for c in EPITOME_COLUMNS])

    cur.execute(f"DELETE FROM [{table_name}]")

    for row in rows:
        values = [row.get(col) for col in EPITOME_COLUMNS]
        cur.execute(
            f"INSERT INTO [{table_name}] ({col_names}) VALUES ({placeholders})",
            values,
        )

    n_success = sum(1 for r in rows if r.get("status") == "success")
    n_error = sum(1 for r in rows if r.get("status") == "error")
    n_skipped = sum(1 for r in rows if r.get("status") == "skipped")

    cur.execute(
        "INSERT INTO run_metadata VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            datetime.now().isoformat(),
            json.dumps(full_config),
            table_name,
            len(rows),
            n_success,
            n_error,
            n_skipped,
        ),
    )

    conn.commit()
    conn.close()
