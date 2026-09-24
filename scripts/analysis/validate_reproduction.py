#!/usr/bin/env python3
"""Validate fixed-cutoff caller and host/virus/sgRNA split outputs."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_BARCODE = {
    "bc04": (1506, 1270, 183),
    "bc05": (11405, 9370, 1201),
    "bc06": (12182, 9799, 1415),
    "bc07": (17284, 14061, 2441),
    "bc11": (44029, 36569, 5722),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--caller-dir", type=Path, required=True)
    parser.add_argument("--split-dir", type=Path, required=True)
    parser.add_argument("--reference-caller", type=Path)
    parser.add_argument("--reference-split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def tables_equal(left_path: Path, right_path: Path) -> bool:
    left = pd.read_csv(left_path, sep="\t")
    right = pd.read_csv(right_path, sep="\t")
    if left.shape != right.shape or list(left.columns) != list(right.columns):
        return False
    for column in left.columns:
        if pd.api.types.is_numeric_dtype(left[column]) and pd.api.types.is_numeric_dtype(right[column]):
            if not np.allclose(
                left[column].to_numpy(float),
                right[column].to_numpy(float),
                rtol=0,
                atol=1e-12,
                equal_nan=True,
            ):
                return False
        elif not left[column].fillna("<NA>").astype(str).equals(
            right[column].fillna("<NA>").astype(str)
        ):
            return False
    return True


def compare_summary_directory(observed: Path, reference: Path) -> int:
    checked = 0
    for reference_path in sorted(reference.glob("*.tsv")):
        observed_path = observed / reference_path.name
        require(observed_path.exists(), f"Missing reproduced table: {observed_path}")
        require(
            tables_equal(observed_path, reference_path),
            f"Summary table differs: {reference_path.name}",
        )
        checked += 1
    return checked


def gzip_contents_equal(left_path: Path, right_path: Path) -> bool:
    with gzip.open(left_path, "rb") as left, gzip.open(right_path, "rb") as right:
        while True:
            left_block = left.read(1024 * 1024)
            right_block = right.read(1024 * 1024)
            if left_block != right_block:
                return False
            if not left_block:
                return True


def main() -> None:
    args = parse_args()
    virus = pd.read_parquet(args.caller_dir / "mhv_reads_terminal_nona_direct.parquet")
    host = pd.read_parquet(args.caller_dir / "host_terminal_nona_direct_scores.parquet")
    require(len(virus) == 86406 and virus["read_id"].nunique() == 86406, "Virus row/ID mismatch")
    require(len(host) == 82709 and host["read_id"].nunique() == 82709, "Host row/ID mismatch")

    eligible = virus["eligible_terminal_nona_direct"].eq(1)
    positive = eligible & virus["terminal_junction_nonA_like"].eq(True)
    require(int(eligible.sum()) == 71069, "Unexpected callable count")
    require(int(positive.sum()) == 10962, "Unexpected non-A-like count")
    require(
        virus.loc[~eligible, "terminal_junction_nonA_like"].isna().all(),
        "Uncallable reads contain binary calls",
    )
    require(
        virus.loc[~eligible, ["terminal_nona_signed_difference_pa", "terminal_nona_abs_difference_pa"]]
        .isna()
        .all()
        .all(),
        "Uncallable reads contain derived scores",
    )

    reference_values = virus["pooled_long_host_A_DNA_reference_pa"].dropna().unique()
    require(len(reference_values) == 1, "Reference is not constant")
    require(np.isclose(reference_values[0], -10.190704345703125, atol=1e-12), "Reference mismatch")
    require(int(host["long_host_reference_member"].sum()) == 26997, "Reference-member count mismatch")

    summary = pd.read_csv(args.caller_dir / "mhv_barcode_nona_direct_summary.tsv", sep="\t")
    for barcode, (total, callable_n, positive_n) in EXPECTED_BARCODE.items():
        row = summary.loc[summary["mapping_barcode"].eq(barcode)]
        require(len(row) == 1, f"Missing barcode row: {barcode}")
        observed = tuple(
            int(row.iloc[0][column])
            for column in ("n_total", "n_callable", "n_terminal_junction_nonA_like")
        )
        require(observed == (total, callable_n, positive_n), f"Barcode count mismatch: {barcode}")

    split_validation = json.loads((args.split_dir / "validation_report.json").read_text())
    expected_split = {
        "n_host_rows": 82709,
        "n_virus_rows": 86406,
        "n_high_confidence_strict_sgrna_reads": 12617,
        "n_high_confidence_strict_sgrna_nonA_like_reads": 1307,
        "n_high_confidence_grna_reads": 52,
        "n_unidentifiable_virus_reads": 73737,
    }
    for key, value in expected_split.items():
        require(split_validation.get(key) == value, f"Split validation mismatch: {key}")
    require(split_validation.get("status") == "PASS", "Split validation did not pass")

    max_core_abs_error = None
    n_reference_summary_tables_checked = 0
    n_reference_gzip_tables_checked = 0
    if args.reference_caller and (args.reference_caller / "mhv_reads_terminal_nona_direct.parquet").exists():
        reference = pd.read_parquet(args.reference_caller / "mhv_reads_terminal_nona_direct.parquet")
        columns = [
            "read_id",
            "eligible_terminal_nona_direct",
            "terminal_nona_direct_qc_status",
            "terminal_minus_stableA_pa",
            "terminal_nona_signed_difference_pa",
            "terminal_nona_abs_difference_pa",
            "terminal_junction_nonA_like",
        ]
        left = virus[columns].sort_values("read_id").reset_index(drop=True)
        right = reference[columns].sort_values("read_id").reset_index(drop=True)
        require(left["read_id"].equals(right["read_id"]), "Reference read IDs differ")
        for column in ("eligible_terminal_nona_direct", "terminal_nona_direct_qc_status", "terminal_junction_nonA_like"):
            require(left[column].equals(right[column]), f"Reference categorical field differs: {column}")
        numeric = ["terminal_minus_stableA_pa", "terminal_nona_signed_difference_pa", "terminal_nona_abs_difference_pa"]
        errors = []
        for column in numeric:
            a = pd.to_numeric(left[column], errors="coerce").to_numpy(float)
            b = pd.to_numeric(right[column], errors="coerce").to_numpy(float)
            require(np.array_equal(np.isnan(a), np.isnan(b)), f"Reference NA pattern differs: {column}")
            finite = np.isfinite(a) & np.isfinite(b)
            errors.append(float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else 0.0)
        max_core_abs_error = max(errors)
        require(max_core_abs_error == 0.0, f"Reference numeric fields differ: {max_core_abs_error}")
        n_reference_summary_tables_checked += compare_summary_directory(
            args.caller_dir, args.reference_caller
        )
        for name in (
            "mhv_reads_terminal_nona_direct.tsv.gz",
            "mhv_terminal_junction_nona_like_reads.tsv.gz",
        ):
            require(
                gzip_contents_equal(args.caller_dir / name, args.reference_caller / name),
                f"Decompressed reference table differs: {name}",
            )
            n_reference_gzip_tables_checked += 1

    if args.reference_split and args.reference_split.exists():
        n_reference_summary_tables_checked += compare_summary_directory(
            args.split_dir, args.reference_split
        )
        name = "virus_high_conf_sgrna_nonA_like_reads.tsv.gz"
        require(
            gzip_contents_equal(args.split_dir / name, args.reference_split / name),
            f"Decompressed reference table differs: {name}",
        )
        n_reference_gzip_tables_checked += 1

    report = {
        "status": "PASS",
        "n_host_rows": len(host),
        "n_virus_rows": len(virus),
        "n_virus_callable": int(eligible.sum()),
        "n_virus_terminal_junction_nonA_like": int(positive.sum()),
        "n_long_host_reference": int(host["long_host_reference_member"].sum()),
        "pooled_reference_pa": float(reference_values[0]),
        "n_strict_sgrna": split_validation["n_high_confidence_strict_sgrna_reads"],
        "n_strict_sgrna_nonA_like": split_validation["n_high_confidence_strict_sgrna_nonA_like_reads"],
        "reference_core_max_abs_error": max_core_abs_error,
        "n_reference_summary_tables_checked": n_reference_summary_tables_checked,
        "n_reference_gzip_tables_checked": n_reference_gzip_tables_checked,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
