#!/usr/bin/env python3
"""Split the fixed-cutoff terminal-junction result by host/virus and MHV sgRNA."""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


BARCODES = ("bc04", "bc05", "bc06", "bc07", "bc11")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-result", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def wilson(successes: int, total: int) -> tuple[float, float]:
    if total == 0:
        return np.nan, np.nan
    z = 1.959963984540054
    p = successes / total
    denominator = 1.0 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def summarise(frame: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    grouper: str | list[str] = group_columns[0] if len(group_columns) == 1 else group_columns
    for keys, part in frame.groupby(grouper, observed=True, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        callable_mask = part["eligible_terminal_nona_direct"].eq(1)
        positive = callable_mask & part["terminal_junction_nonA_like"].eq(True)
        n_total = len(part)
        n_callable = int(callable_mask.sum())
        n_positive = int(positive.sum())
        low, high = wilson(n_positive, n_callable)
        row = {column: value for column, value in zip(group_columns, keys, strict=True)}
        row.update(
            {
                "n_total": n_total,
                "n_callable": n_callable,
                "callable_fraction": n_callable / n_total if n_total else np.nan,
                "n_terminal_junction_nonA_like": n_positive,
                "n_A_junction_like": n_callable - n_positive,
                "nonA_like_fraction_callable": n_positive / n_callable if n_callable else np.nan,
                "wilson95_low_descriptive": low,
                "wilson95_high_descriptive": high,
                "all_read_lower_bound": n_positive / n_total if n_total else np.nan,
                "median_abs_difference_pa_callable": part.loc[
                    callable_mask, "terminal_nona_abs_difference_pa"
                ].median(),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def length_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    callable_frame = frame.loc[
        frame["sgrna_primary_high_confidence"].eq(1)
        & frame["eligible_terminal_nona_direct"].eq(1)
    ].copy()
    for keys, part in callable_frame.groupby(
        [
            "sample_group",
            "mapping_barcode",
            "sgrna_primary_assignment",
            "terminal_nona_direct_class",
        ],
        observed=True,
        dropna=False,
        sort=True,
    ):
        row: dict[str, object] = {
            "sample_group": keys[0],
            "mapping_barcode": keys[1],
            "sgrna_primary_assignment": keys[2],
            "terminal_nona_direct_class": keys[3],
            "n_reads": len(part),
        }
        for source, prefix in (("polya_length_nt", "polya_nt"), ("read_nt_len", "read_nt")):
            values = pd.to_numeric(part[source], errors="coerce").dropna()
            row[f"{prefix}_n"] = len(values)
            row[f"{prefix}_median"] = values.median() if len(values) else np.nan
            row[f"{prefix}_q25"] = values.quantile(0.25) if len(values) else np.nan
            row[f"{prefix}_q75"] = values.quantile(0.75) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    host_path = args.source_result / "host_terminal_nona_direct_scores.parquet"
    virus_path = args.source_result / "mhv_reads_terminal_nona_direct.parquet"
    host = pd.read_parquet(host_path)
    virus = pd.read_parquet(virus_path)

    if len(host) != 82_709 or len(virus) != 86_406:
        raise RuntimeError(f"Unexpected source row counts: host={len(host)}, virus={len(virus)}")
    if host["read_id"].duplicated().any() or virus["read_id"].duplicated().any():
        raise RuntimeError("Duplicate read IDs in source results")

    host = host.copy()
    virus = virus.copy()
    host["read_origin"] = "host_mapped"
    virus["read_origin"] = "MHV_mapped"
    combined = pd.concat([host, virus[host.columns]], ignore_index=True)

    barcode_origin = summarise(
        combined,
        ["sample_group", "mapping_barcode", "read_origin"],
    )
    barcode_origin["pooled_reference_pa"] = float(
        combined["pooled_long_host_A_DNA_reference_pa"].dropna().iloc[0]
    )
    barcode_origin["fixed_cutoff_pa"] = float(
        combined["terminal_nona_fixed_cutoff_pa"].dropna().iloc[0]
    )

    high_conf = virus.loc[virus["sgrna_primary_high_confidence"].eq(1)].copy()
    grna = high_conf.loc[high_conf["sgrna_primary_assignment"].eq("gRNA_ORF1ab")].copy()
    strict_sgrna = high_conf.loc[
        ~high_conf["sgrna_primary_assignment"].eq("gRNA_ORF1ab")
    ].copy()
    unassigned = virus.loc[
        virus["sgrna_all_read_identifiability"].eq("nested_RNA_unidentifiable")
    ].copy()
    sgrna_barcode = summarise(
        strict_sgrna,
        ["sample_group", "mapping_barcode", "sgrna_primary_assignment"],
    )
    sgrna_group = summarise(
        strict_sgrna,
        ["sample_group", "sgrna_primary_assignment"],
    )
    grna_barcode = summarise(
        grna,
        ["sample_group", "mapping_barcode", "sgrna_primary_assignment"],
    )
    unassigned_barcode = summarise(
        unassigned,
        ["sample_group", "mapping_barcode", "sgrna_all_read_identifiability"],
    )
    identifiability = summarise(
        virus,
        ["sample_group", "mapping_barcode", "sgrna_all_read_identifiability"],
    )
    lengths = length_summary(strict_sgrna)

    positive_high_conf = strict_sgrna.loc[
        strict_sgrna["eligible_terminal_nona_direct"].eq(1)
        & strict_sgrna["terminal_junction_nonA_like"].eq(True)
    ].copy()
    positive_columns = [
        "read_id",
        "sample_group",
        "mapping_barcode",
        "sgrna_primary_assignment",
        "sgrna_assignment_method",
        "mapping_mapq",
        "dorado_mean_qscore",
        "terminal_nona_signed_difference_pa",
        "terminal_nona_abs_difference_pa",
        "terminal_nona_fixed_cutoff_pa",
        "polya_length_nt",
        "polya_length_source",
        "read_nt_len",
    ]
    positive_high_conf = positive_high_conf[positive_columns].sort_values(
        ["sample_group", "mapping_barcode", "sgrna_primary_assignment", "read_id"]
    )

    write_tsv(barcode_origin, args.outdir / "barcode_host_vs_virus_nona_summary.tsv")
    write_tsv(sgrna_barcode, args.outdir / "virus_high_conf_sgrna_by_barcode_summary.tsv")
    write_tsv(sgrna_group, args.outdir / "virus_high_conf_sgrna_by_group_summary.tsv")
    write_tsv(grna_barcode, args.outdir / "virus_high_conf_grna_by_barcode_summary.tsv")
    write_tsv(unassigned_barcode, args.outdir / "virus_unassigned_by_barcode_summary.tsv")
    write_tsv(identifiability, args.outdir / "virus_transcript_identifiability_by_barcode.tsv")
    write_tsv(lengths, args.outdir / "virus_high_conf_sgrna_length_by_barcode.tsv")
    positive_high_conf.to_csv(
        args.outdir / "virus_high_conf_sgrna_nonA_like_reads.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    accounting_rows: list[dict[str, object]] = []
    for barcode in BARCODES:
        all_part = virus.loc[virus["mapping_barcode"].eq(barcode)]
        unassigned_part = unassigned.loc[unassigned["mapping_barcode"].eq(barcode)]
        grna_part = grna.loc[grna["mapping_barcode"].eq(barcode)]
        sgrna_part = strict_sgrna.loc[strict_sgrna["mapping_barcode"].eq(barcode)]
        row: dict[str, object] = {"mapping_barcode": barcode}
        for label, part in (
            ("all", all_part),
            ("unassigned", unassigned_part),
            ("gRNA", grna_part),
            ("strict_sgRNA", sgrna_part),
        ):
            callable_mask = part["eligible_terminal_nona_direct"].eq(1)
            positive_mask = callable_mask & part["terminal_junction_nonA_like"].eq(True)
            row[f"{label}_total"] = len(part)
            row[f"{label}_callable"] = int(callable_mask.sum())
            row[f"{label}_nonA_like"] = int(positive_mask.sum())
        row["total_reconciles"] = row["all_total"] == (
            row["unassigned_total"] + row["gRNA_total"] + row["strict_sgRNA_total"]
        )
        row["callable_reconciles"] = row["all_callable"] == (
            row["unassigned_callable"] + row["gRNA_callable"] + row["strict_sgRNA_callable"]
        )
        row["nonA_like_reconciles"] = row["all_nonA_like"] == (
            row["unassigned_nonA_like"] + row["gRNA_nonA_like"] + row["strict_sgRNA_nonA_like"]
        )
        accounting_rows.append(row)
    accounting = pd.DataFrame(accounting_rows)
    write_tsv(accounting, args.outdir / "virus_transcript_accounting_by_barcode.tsv")

    validation = {
        "status": "PASS",
        "n_host_rows": len(host),
        "n_virus_rows": len(virus),
        "n_barcode_origin_summary_total": int(barcode_origin["n_total"].sum()),
        "n_high_confidence_transcript_reads": len(high_conf),
        "n_high_confidence_strict_sgrna_reads": len(strict_sgrna),
        "n_high_confidence_grna_reads": len(grna),
        "n_unidentifiable_virus_reads": int(
            virus["sgrna_all_read_identifiability"].eq("nested_RNA_unidentifiable").sum()
        ),
        "n_high_confidence_strict_sgrna_nonA_like_reads": len(positive_high_conf),
        "n_positive_reads_missing_polya_length": int(positive_high_conf["polya_length_nt"].isna().sum()),
        "n_positive_reads_missing_read_length": int(positive_high_conf["read_nt_len"].isna().sum()),
        "high_conf_sgrna_barcode_total_matches": int(sgrna_barcode["n_total"].sum()) == len(strict_sgrna),
        "high_conf_sgrna_positive_matches": int(
            sgrna_barcode["n_terminal_junction_nonA_like"].sum()
        ) == len(positive_high_conf),
        "identifiability_total_matches": int(identifiability["n_total"].sum()) == len(virus),
        "transcript_accounting_all_barcodes_pass": bool(
            accounting[
                ["total_reconciles", "callable_reconciles", "nonA_like_reconciles"]
            ].to_numpy().all()
        ),
    }
    if not all(
        validation[key]
        for key in (
            "high_conf_sgrna_barcode_total_matches",
            "high_conf_sgrna_positive_matches",
            "identifiability_total_matches",
            "transcript_accounting_all_barcodes_pass",
        )
    ):
        raise RuntimeError(f"Validation failed: {validation}")
    (args.outdir / "validation_report.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )

    readme = """# Host/virus and MHV sgRNA split of the fixed-cutoff non-A-like result

`barcode_host_vs_virus_nona_summary.tsv` separates host-mapped and MHV-mapped
reads within every barcode. Host reads are not assigned to MHV sgRNAs.

The two `virus_high_conf_sgrna_*` summary files use only diagnostic sgRNA
assignments supported by canonical leader-body junctions. The conservative
ORF1 gRNA calls are separate in `virus_high_conf_grna_by_barcode_summary.tsv`.
`virus_transcript_identifiability_by_barcode.tsv` retains all
MHV-mapped reads and explicitly labels non-diagnostic short/nested reads as
`nested_RNA_unidentifiable`; these reads are not forced into N or another gene.

`virus_high_conf_sgrna_nonA_like_reads.tsv.gz` contains strict high-confidence
sgRNA reads exceeding the fixed 5 pA cutoff, including poly(A) length and read
length. `virus_unassigned_by_barcode_summary.tsv` and
`virus_transcript_accounting_by_barcode.tsv` make the unassigned and total-read
accounting explicit. The endpoint remains an operational non-A-like junction deviation, not
a standards-validated biochemical U call.
"""
    (args.outdir / "README.md").write_text(readme)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "python": platform.python_version(),
        "script": str(Path(__file__).resolve()),
        "source_result": str(args.source_result.resolve()),
    }
    (args.outdir / "result_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
