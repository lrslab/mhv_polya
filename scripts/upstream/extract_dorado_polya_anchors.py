#!/usr/bin/env python3
"""Extract Dorado poly(A) anchor coordinates from unaligned BAM files.

The Dorado ``pa:B:i`` tag stores, in signal-sample coordinates, the search
anchor, primary poly(A) start/end, and optional secondary start/end.  Keeping
these values in a small parquet table avoids repeatedly scanning the large BAM
files and lets downstream callers join them to WarpDemuX-aligned POD5 windows.
"""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


BARCODES = ("bc04", "bc05", "bc06", "bc07", "bc11")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam-dir",
        type=Path,
        required=True,
        help="Directory containing bc04.polya.bam, ..., bc11.polya.bam",
    )
    parser.add_argument(
        "--boundary-dir",
        type=Path,
        help="Optional WarpDemuX directory containing detected_boundaries_*.csv",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def extract_one(path: Path, barcode: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        for record in bam.fetch(until_eof=True):
            pa = list(record.get_tag("pa")) if record.has_tag("pa") else []
            pa = (pa + [np.nan] * 5)[:5]
            rows.append(
                {
                    "read_id": record.query_name,
                    "dorado_barcode": barcode,
                    "dorado_polya_estimate_nt": (
                        record.get_tag("pt") if record.has_tag("pt") else np.nan
                    ),
                    "dorado_pa_anchor": pa[0],
                    "dorado_polya_start": pa[1],
                    "dorado_polya_end": pa[2],
                    "dorado_secondary_polya_start": pa[3],
                    "dorado_secondary_polya_end": pa[4],
                    "dorado_trimmed_signal_start": (
                        record.get_tag("ts") if record.has_tag("ts") else np.nan
                    ),
                    "dorado_basecalled_signal_end": (
                        record.get_tag("ns") if record.has_tag("ns") else np.nan
                    ),
                    "dorado_scaling_midpoint_pa": (
                        record.get_tag("sm") if record.has_tag("sm") else np.nan
                    ),
                    "dorado_scaling_dispersion_pa": (
                        record.get_tag("sd") if record.has_tag("sd") else np.nan
                    ),
                    "dorado_mean_qscore": (
                        record.get_tag("qs") if record.has_tag("qs") else np.nan
                    ),
                    "dorado_source_bam": path.name,
                }
            )
    frame = pd.DataFrame(rows)
    integer_columns = [
        "dorado_polya_estimate_nt",
        "dorado_pa_anchor",
        "dorado_polya_start",
        "dorado_polya_end",
        "dorado_secondary_polya_start",
        "dorado_secondary_polya_end",
        "dorado_trimmed_signal_start",
        "dorado_basecalled_signal_end",
    ]
    for column in integer_columns:
        frame[column] = pd.array(frame[column], dtype="Int64")
    return frame


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    boundary = None
    if args.boundary_dir is not None:
        boundary_files = sorted(
            glob.glob(str(args.boundary_dir / "detected_boundaries_[0-9]*.csv"))
        )
        if not boundary_files:
            raise FileNotFoundError(f"No boundary CSV files in {args.boundary_dir}")
        boundary_columns = [
            "read_id",
            "adapter_end",
            "polya_start",
            "polya_end",
            "adapter_dt_med",
            "adapter_dt_mad",
            "cnn_adapter_end",
            "llr_adapter_end",
        ]
        boundary = pd.concat(
            [pd.read_csv(path, usecols=boundary_columns) for path in boundary_files],
            ignore_index=True,
        )
        if boundary["read_id"].duplicated().any():
            raise RuntimeError("WarpDemuX boundary read IDs are not unique")
        boundary = boundary.rename(
            columns={
                "adapter_end": "warpdemux_adapter_end",
                "polya_start": "warpdemux_polya_start",
                "polya_end": "warpdemux_polya_end",
                "adapter_dt_med": "warpdemux_adapter_dwell_median_samples",
                "adapter_dt_mad": "warpdemux_adapter_dwell_mad_samples",
                "cnn_adapter_end": "warpdemux_cnn_adapter_end",
                "llr_adapter_end": "warpdemux_llr_adapter_end",
            }
        )
    summary = []
    all_ids: set[str] = set()
    for barcode in BARCODES:
        bam_path = args.bam_dir / f"{barcode}.polya.bam"
        if not bam_path.exists():
            raise FileNotFoundError(bam_path)
        frame = extract_one(bam_path, barcode)
        if boundary is not None:
            frame = frame.merge(boundary, on="read_id", how="left", validate="one_to_one")
            frame["warpdemux_boundary_method"] = np.where(
                frame["warpdemux_cnn_adapter_end"].notna(),
                "cnn",
                np.where(frame["warpdemux_llr_adapter_end"].notna(), "llr", "missing"),
            )
        duplicated = int(frame["read_id"].duplicated().sum())
        if duplicated:
            raise RuntimeError(f"{bam_path} contains {duplicated} duplicate read IDs")
        overlap = all_ids.intersection(frame["read_id"])
        if overlap:
            raise RuntimeError(f"Read IDs occur in multiple barcode BAMs: {len(overlap)}")
        all_ids.update(frame["read_id"])
        frame.to_parquet(args.outdir / f"{barcode}.dorado_polya_anchors.parquet", index=False)
        valid = (
            frame["dorado_polya_estimate_nt"].gt(0)
            & frame["dorado_polya_start"].ge(0)
            & frame["dorado_polya_end"].ge(0)
        )
        summary.append(
            {
                "barcode": barcode,
                "n_records": len(frame),
                "n_valid_primary_polya": int(valid.sum()),
                "n_anchor_equals_ts": int(
                    frame["dorado_pa_anchor"]
                    .eq(frame["dorado_trimmed_signal_start"])
                    .sum()
                ),
            }
        )
    pd.DataFrame(summary).to_csv(
        args.outdir / "dorado_polya_anchor_summary.tsv", sep="\t", index=False
    )
    (args.outdir / "source_definition.json").write_text(
        json.dumps(
            {
                "pa_tag_order": [
                    "polyA_search_anchor",
                    "primary_polyA_start",
                    "primary_polyA_end",
                    "secondary_polyA_start",
                    "secondary_polyA_end",
                ],
                "coordinate_unit": "raw signal samples",
                "barcodes": list(BARCODES),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
