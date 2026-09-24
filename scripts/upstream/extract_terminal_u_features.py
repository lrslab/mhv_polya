#!/usr/bin/env python3
"""Extract calibrated junction and internal poly(A) signals for terminal-U calling.

The source project is read-only. This script reads WarpDemux boundaries and the
barcode POD5 files, converts signal to pA through ``ReadRecord.signal_pa``, and
writes compact fixed-width feature matrices under a user-writable output root.
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pod5


BARCODES = ("bc04", "bc05", "bc06", "bc07", "bc11")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/t3/qqjiang_x112/MHV_polyA"),
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--window", type=int, default=360)
    parser.add_argument("--internal-offset", type=int, default=300)
    parser.add_argument("--baseline-end", type=int, default=900)
    return parser.parse_args()


def load_metadata(project_root: Path) -> pd.DataFrame:
    path = (
        project_root
        / "polyA_host_virus_by_barcode_warpdemux"
        / "all_reads.warpdemux_host_virus_other.polya.tsv"
    )
    usecols = [
        "read_id",
        "warpdemux_barcode",
        "warpdemux_barcode_confidence",
        "mapping_group",
        "mapping_barcode",
        "barcode_mapping_consistent",
        "mapping_status_detail",
        "mapping_reference",
        "mapping_mapq",
        "warpdemux_polya_len_samples",
        "read_nt_len",
        "warpdemux_polya_estimate_nt",
        "warpdemux_signal_len",
        "boundary_source_file",
        "warpdemux_status",
    ]
    frame = pd.read_csv(path, sep="\t", usecols=usecols)
    frame = frame[
        frame["mapping_group"].isin(["host", "virus"])
        & frame["mapping_barcode"].isin(BARCODES)
        & frame["warpdemux_status"].eq("success")
    ].copy()
    if frame["read_id"].duplicated().any():
        raise RuntimeError("Duplicate read IDs in host/virus metadata")
    return frame


def load_boundaries(project_root: Path, wanted: set[str]) -> pd.DataFrame:
    boundary_dir = (
        project_root
        / "MHV_Multiple_sample"
        / "Warpdemux_polyA"
        / "warpdemux_WDX6_rna004_v1_0_20260318_1039_ade68d6a"
        / "boundaries"
    )
    usecols = [
        "read_id",
        "signal_len",
        "adapter_start",
        "adapter_end",
        "polya_start",
        "polya_end",
        "polya_len",
        "polya_mean",
        "polya_std",
        "polya_med",
        "polya_mad",
        "llr_adapter_end",
        "llr_polya_end",
        "cnn_adapter_end",
        "cnn_polya_end",
        "start_peak_adapter_end",
        "start_peak_polya_end",
    ]
    chunks: list[pd.DataFrame] = []
    for path in sorted(boundary_dir.glob("detected_boundaries_*.csv")):
        part = pd.read_csv(path, usecols=usecols)
        part = part[part["read_id"].isin(wanted)].copy()
        part["boundary_file"] = path.name
        chunks.append(part)
    if not chunks:
        raise RuntimeError(f"No boundary files found under {boundary_dir}")
    frame = pd.concat(chunks, ignore_index=True)
    if frame["read_id"].duplicated().any():
        raise RuntimeError("Duplicate read IDs in WarpDemux boundaries")
    return frame


def robust_sigma(values: np.ndarray) -> float:
    med = float(np.nanmedian(values))
    mad = float(np.nanmedian(np.abs(values - med)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma < 0.25:
        sigma = float(np.nanstd(values))
    return sigma


def fixed_slice(signal: np.ndarray, start: int, width: int) -> np.ndarray:
    result = np.full(width, np.nan, dtype=np.float32)
    lo = max(0, int(start))
    hi = min(len(signal), lo + width)
    if hi > lo:
        result[: hi - lo] = signal[lo:hi]
    return result


def extract_barcode(
    frame: pd.DataFrame,
    pod5_path: Path,
    outdir: Path,
    window: int,
    internal_offset: int,
    baseline_end: int,
) -> None:
    frame = frame.reset_index(drop=True).copy()
    n_reads = len(frame)
    junction = np.full((n_reads, window), np.nan, dtype=np.float32)
    internal = np.full((n_reads, window), np.nan, dtype=np.float32)
    lookup = {read_id: i for i, read_id in enumerate(frame["read_id"])}

    baseline_pa = np.full(n_reads, np.nan, dtype=np.float32)
    baseline_sigma_pa = np.full(n_reads, np.nan, dtype=np.float32)
    actual_polya_len = np.zeros(n_reads, dtype=np.int32)
    signal_found = np.zeros(n_reads, dtype=np.int8)
    boundary_valid = np.zeros(n_reads, dtype=np.int8)

    logging.info("Opening %s for %d selected reads", pod5_path, n_reads)
    with pod5.Reader(pod5_path) as reader:
        for count, record in enumerate(reader.reads(selection=list(lookup)), start=1):
            read_id = str(record.read_id)
            idx = lookup.get(read_id)
            if idx is None:
                continue
            signal_found[idx] = 1
            row = frame.iloc[idx]
            signal = np.asarray(record.signal_pa, dtype=np.float32)
            start = int(row["polya_start"])
            end = min(int(row["polya_end"]), len(signal))
            actual_polya_len[idx] = max(0, end - start)
            if start < 0 or start >= len(signal) or end <= start:
                continue

            base_start = start + internal_offset
            base_stop = min(start + baseline_end, end - 60)
            if base_stop - base_start < 180:
                continue
            base = signal[base_start:base_stop]
            baseline = float(np.nanmedian(base))
            sigma = robust_sigma(base)
            if not np.isfinite(baseline) or not np.isfinite(sigma) or sigma <= 0:
                continue

            baseline_pa[idx] = baseline
            baseline_sigma_pa[idx] = sigma
            junction[idx] = fixed_slice(signal, start, window) - baseline
            internal[idx] = fixed_slice(signal, base_start, window) - baseline
            boundary_valid[idx] = 1

            if count % 10_000 == 0:
                logging.info("%s: read %d", pod5_path.name, count)

    frame["signal_found"] = signal_found
    frame["boundary_valid"] = boundary_valid
    frame["actual_polya_len_samples"] = actual_polya_len
    frame["baseline_pa"] = baseline_pa
    frame["baseline_sigma_pa"] = baseline_sigma_pa
    frame["callable_base"] = (
        frame["boundary_valid"].eq(1)
        & frame["actual_polya_len_samples"].ge(baseline_end)
        & np.isfinite(frame["baseline_pa"])
        & np.isfinite(frame["baseline_sigma_pa"])
    ).astype(np.int8)
    frame["boundary_method"] = np.where(
        frame["cnn_adapter_end"].notna() & frame["cnn_polya_end"].notna(),
        "cnn",
        "fallback",
    )

    barcode = str(frame["mapping_barcode"].iloc[0])
    np.save(outdir / f"{barcode}.junction_pa_centered.npy", junction)
    np.save(outdir / f"{barcode}.internal_pa_centered.npy", internal)
    frame.to_parquet(outdir / f"{barcode}.metadata.parquet", index=False)
    logging.info(
        "%s complete: found=%d/%d valid=%d callable=%d",
        barcode,
        int(signal_found.sum()),
        n_reads,
        int(boundary_valid.sum()),
        int(frame["callable_base"].sum()),
    )


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    os.chdir(args.outdir)  # POD5 v3/v4 migration needs a writable temp location.
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    metadata = load_metadata(args.project_root)
    boundaries = load_boundaries(args.project_root, set(metadata["read_id"]))
    merged = metadata.merge(boundaries, on="read_id", how="left", validate="one_to_one")
    if merged["polya_start"].isna().any():
        missing = int(merged["polya_start"].isna().sum())
        raise RuntimeError(f"{missing} selected reads lack a WarpDemux boundary")

    pod5_dir = (
        args.project_root
        / "MHV_Multiple_sample"
        / "20260306_5bar"
        / "20260306_1823_MN25294_FAZ51928_089bb3a7"
        / "pod5_by_barcode"
    )
    for barcode in BARCODES:
        subset = merged[merged["mapping_barcode"].eq(barcode)].copy()
        extract_barcode(
            subset,
            pod5_dir / f"{barcode}.pod5",
            args.outdir,
            args.window,
            args.internal_offset,
            args.baseline_end,
        )


if __name__ == "__main__":
    main()
