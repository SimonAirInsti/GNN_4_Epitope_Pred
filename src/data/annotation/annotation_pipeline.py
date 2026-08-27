#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
annotation_pipeline.py — SAbDab2 epitope annotation pipeline.

Reads CIF files for antibody-antigen complexes, computes contact distances
and RSA, labels epitope residues, and writes results to a SQLite database.

Usage:
    cd <PROJECT_ROOT>
    python src/data/annotation/annotation_pipeline.py
"""

import os
import sys
import logging
import json
from datetime import datetime
from pathlib import Path

# Ensure project root is on sys.path so 'from src...' imports work regardless
# of the working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data.annotation.config_loader import load_annotation_config
from src.data.annotation.csv_loader import load_csv
from src.data.annotation.cif_parser import load_complex_coords
from src.data.annotation.rsa_calculator import compute_rsa
from src.data.annotation.contact_calculator import compute_distances
from src.data.annotation.epitope_labeler import (
    label_epitope,
    format_epitope_labels,
    format_paratope_labels,
    count_epitope_labels,
    count_paratope_labels,
    serialize_sequences,
)
from src.data.annotation.db_writer import create_database, insert_rows


def _setup_logging(log_cfg, log_dir):
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"annotation_{timestamp}.log")

    formatter = logging.Formatter(
        "[%(levelname)s] %(asctime)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)

    handlers = [console]
    if log_cfg.get("log_to_file", True):
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)

    logger = logging.getLogger("annotation")
    logger.setLevel(getattr(logging, log_cfg.get("level", "INFO")))
    for h in handlers:
        logger.addHandler(h)

    return logger


def _build_row(row, result, epitope_labels, paratope_dists,
               contact_pairs, ann_cfg):
    """Build a single DB row dict from annotation results."""
    cutoff = ann_cfg["contact_cutoff"]
    rsa_thresh = ann_cfg["rsa_threshold"]

    epitope_str = format_epitope_labels(epitope_labels)
    paratope_str = format_paratope_labels(paratope_dists, cutoff=cutoff)
    n_epi = count_epitope_labels(epitope_labels)
    n_para = count_paratope_labels(paratope_dists, cutoff=cutoff)
    n_ag = sum(
        len(res_dict)
        for res_dict in result["antigen_coords"].values()
    )

    return {
        "INSTANCE":             row["INSTANCE"],
        "PDB_ID":               row.get("PDB_ID"),
        "SABDAB_ID":            row.get("SABDAB_ID"),
        "HEAVY_ID":             row.get("HEAVY_ID"),
        "LIGHT_ID":             row.get("LIGHT_ID"),
        "Hchain":               row.get("Hchain"),
        "Lchain":               row.get("Lchain"),
        "agchains":             row.get("agchains"),
        "agtypes":              row.get("agtypes"),
        "agresolvedseqs":       row.get("agresolvedseqs"),
        "agexpectedseqs":       row.get("agexpectedseqs"),
        "cdrh3_cluster":        row.get("cdrh3_cluster"),
        "cdrh123_cluster":      row.get("cdrh123_cluster"),
        "cdrl123_cluster":      row.get("cdrl123_cluster"),
        "ab_cluster":           row.get("ab_cluster"),
        "agclusters":           row.get("agclusters"),
        "ab_ag_cluster":        row.get("ab_ag_cluster"),
        "ab_ag_split":          row.get("ab_ag_split"),
        "antigen_sequences":    serialize_sequences(result["antigen_sequences"]),
        "antibody_sequences":   serialize_sequences(result["antibody_sequences"]),
        "epitope_labels":       epitope_str,
        "paratope_labels":      paratope_str,
        "n_epitope_residues":   n_epi,
        "n_paratope_residues":  n_para,
        "n_contact_pairs":      len(contact_pairs),
        "n_ag_residues":        n_ag,
        "contact_cutoff":       cutoff,
        "rsa_threshold":        rsa_thresh,
        "status":               "success",
        "error_message":        None,
    }


def run(config_path=None):
    """Main pipeline entry point."""
    cfg = load_annotation_config(config_path)

    project_root = Path(__file__).resolve().parents[3]
    csv_path = project_root / cfg["paths"]["csv_path"]
    cif_dir = project_root / cfg["paths"]["cif_dir"]
    db_dir = project_root / cfg["paths"]["output_db_dir"]
    db_path = db_dir / "db_sabdab_epitopes.db"

    os.makedirs(db_dir, exist_ok=True)

    logger = _setup_logging(cfg["logging"], str(db_dir))
    logger.info("Annotation pipeline started")
    logger.info(f"Config loaded: cutoff={cfg['annotation']['contact_cutoff']}, "
                f"rsa_threshold={cfg['annotation']['rsa_threshold']}, "
                f"min_chain_residues={cfg['annotation']['min_chain_residues']}")

    # --- Load + Filter CSV ---
    logger.info(f"Loading CSV: {csv_path}")
    raw_df = load_csv(str(csv_path))
    logger.info(f"CSV loaded: {len(raw_df)} rows")

    logger.info("Applying filters...")
    df = load_csv(str(csv_path), filter_cfg=cfg["filter"])
    logger.info(f"Filter summary: {len(raw_df)} -> {len(df)} complexes")
    logger.info(f"Complexes to annotate: {len(df)}")

    # --- Create DB ---
    ann_cfg = cfg["annotation"]
    table_name = create_database(str(db_path), ann_cfg)
    logger.info(f"Database created: {db_path}, table: {table_name}")

    # --- Annotate ---
    interval = cfg["logging"].get("progress_interval", 50)
    rows = []
    success_count = 0
    error_count = 0
    skip_count = 0

    total = len(df)
    for i, (_, row) in enumerate(df.iterrows()):
        instance = row["INSTANCE"]
        if (i + 1) % interval == 0 or i == 0:
            pct = (i + 1) / total * 100
            avg_epi = (
                sum(r["n_epitope_residues"] for r in rows if r["status"] == "success")
                / max(success_count, 1)
            )
            logger.info(
                f"Progress: {i + 1}/{total} ({pct:.1f}%) "
                f"— avg epitope residues: {avg_epi:.1f}"
            )

        try:
            result = load_complex_coords(
                row, str(cif_dir),
                min_chain_residues=ann_cfg["min_chain_residues"],
                logger=logger,
            )

            rsa_all = compute_rsa(
                result["cif_path"],
                point_number=ann_cfg["sasa_point_number"],
            )
            rsa_ag = {
                ch: rsa_all.get(ch, {})
                for ch in result["antigen_chains"]
            }

            ag_to_ab = compute_distances(
                result["antigen_coords"],
                result["antibody_coords"],
                cutoff=ann_cfg["contact_cutoff"],
            )
            ab_to_ag = compute_distances(
                result["antibody_coords"],
                result["antigen_coords"],
                cutoff=ann_cfg["contact_cutoff"],
            )

            epitope_labels = label_epitope(
                result["antigen_coords"],
                ag_to_ab["query_dists"],
                rsa_ag,
                cutoff=ann_cfg["contact_cutoff"],
                rsa_threshold=ann_cfg["rsa_threshold"],
            )

            db_row = _build_row(
                row, result, epitope_labels,
                ab_to_ag["query_dists"],
                ag_to_ab["contact_pairs"],
                ann_cfg,
            )
            rows.append(db_row)
            success_count += 1

        except Exception as e:
            error_count += 1
            logger.error(f"Failed {instance}: {e}")
            rows.append({
                "INSTANCE": instance,
                "PDB_ID": row.get("PDB_ID"),
                "SABDAB_ID": row.get("SABDAB_ID"),
                "HEAVY_ID": row.get("HEAVY_ID"),
                "LIGHT_ID": row.get("LIGHT_ID"),
                "Hchain": row.get("Hchain"),
                "Lchain": row.get("Lchain"),
                "agchains": row.get("agchains"),
                "agtypes": row.get("agtypes"),
                "agresolvedseqs": row.get("agresolvedseqs"),
                "agexpectedseqs": row.get("agexpectedseqs"),
                "cdrh3_cluster": row.get("cdrh3_cluster"),
                "cdrh123_cluster": row.get("cdrh123_cluster"),
                "cdrl123_cluster": row.get("cdrl123_cluster"),
                "ab_cluster": row.get("ab_cluster"),
                "agclusters": row.get("agclusters"),
                "ab_ag_cluster": row.get("ab_ag_cluster"),
                "ab_ag_split": row.get("ab_ag_split"),
                "antigen_sequences": None,
                "antibody_sequences": None,
                "epitope_labels": None,
                "paratope_labels": None,
                "n_epitope_residues": 0,
                "n_paratope_residues": 0,
                "n_contact_pairs": 0,
                "n_ag_residues": 0,
                "contact_cutoff": ann_cfg["contact_cutoff"],
                "rsa_threshold": ann_cfg["rsa_threshold"],
                "status": "error",
                "error_message": str(e),
            })

    # --- Write DB ---
    logger.info(f"Writing {len(rows)} rows to database...")
    insert_rows(str(db_path), table_name, rows, cfg)

    # --- Summary ---
    success_rows = [r for r in rows if r["status"] == "success"]
    epi_counts = [r["n_epitope_residues"] for r in success_rows]

    logger.info("=" * 60)
    logger.info("ANNOTATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total processed:  {len(rows)}")
    logger.info(f"  Success:        {success_count}")
    logger.info(f"  Error:          {error_count}")
    logger.info(f"  Skipped:        {skip_count}")
    logger.info(f"DB path:          {db_path}")
    logger.info(f"Table name:       {table_name}")

    if epi_counts:
        logger.info(f"Epitope residues per complex:")
        logger.info(f"  Mean:    {sum(epi_counts) / len(epi_counts):.1f}")
        epi_sorted = sorted(epi_counts)
        median = epi_sorted[len(epi_sorted) // 2]
        logger.info(f"  Median:  {median}")
        logger.info(f"  Min:     {min(epi_counts)}")
        logger.info(f"  Max:     {max(epi_counts)}")
        zero_epi = sum(1 for c in epi_counts if c == 0)
        logger.info(f"  Zero epitopes: {zero_epi} / {len(epi_counts)}")

    if error_count > 0:
        error_instances = [
            r["INSTANCE"] for r in rows if r["status"] == "error"
        ]
        logger.info(f"Failed instances ({error_count}): {error_instances[:20]}")
        if len(error_instances) > 20:
            logger.info(f"  ... and {len(error_instances) - 20} more")

    logger.info("Annotation pipeline completed")
    return str(db_path), table_name


if __name__ == "__main__":
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else None
    db, table = run(cfg_path)
    print(f"\nDatabase: {db}")
    print(f"Table:    {table}")
