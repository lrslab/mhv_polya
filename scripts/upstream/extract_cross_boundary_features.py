#!/usr/bin/env python3
"""Extract calibrated signal on both sides of WarpDemuX polyA_start.

For a sequence 5'-...A^L U^n-3', direct-RNA signal arrives 3'->5'.  The
terminal U states therefore occur on the adapter-facing side of the A plateau.
This extractor retains samples before and after the boundary so they cannot be
silently clipped as in a one-sided poly(A) slice.
"""

from __future__ import annotations

import argparse
import json
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
        "--project-root", type=Path, default=Path("/t3/qqjiang_x112/MHV_polyA")
    )
    parser.add_argument("--base-features", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pre-boundary", type=int, default=240)
    parser.add_argument("--post-boundary", type=int, default=360)
    return parser.parse_args()


def fixed_slice(signal: np.ndarray, start: int, width: int) -> np.ndarray:
    result = np.full(width, np.nan, dtype=np.float32)
    source_lo = max(0, start)
    source_hi = min(len(signal), start + width)
    if source_hi <= source_lo:
        return result
    destination_lo = source_lo - start
    destination_hi = destination_lo + source_hi - source_lo
    result[destination_lo:destination_hi] = signal[source_lo:source_hi]
    return result


def extract_barcode(
    barcode: str,
    metadata: pd.DataFrame,
    pod5_path: Path,
    outdir: Path,
    pre_boundary: int,
    post_boundary: int,
) -> None:
    metadata = metadata.reset_index(drop=True).copy()
    width = pre_boundary + post_boundary
    matrix = np.full((len(metadata), width), np.nan, dtype=np.float32)
    found = np.zeros(len(metadata), dtype=np.int8)
    complete_window = np.zeros(len(metadata), dtype=np.int8)
    lookup = {read_id: idx for idx, read_id in enumerate(metadata["read_id"])}

    with pod5.Reader(pod5_path) as reader:
        for count, record in enumerate(reader.reads(selection=list(lookup)), start=1):
            idx = lookup.get(str(record.read_id))
            if idx is None:
                continue
            found[idx] = 1
            baseline = float(metadata.iloc[idx]["baseline_pa"])
            if not np.isfinite(baseline):
                continue
            boundary = int(metadata.iloc[idx]["polya_start"])
            signal = np.asarray(record.signal_pa, dtype=np.float32)
            start = boundary - pre_boundary
            matrix[idx] = fixed_slice(signal, start, width) - baseline
            complete_window[idx] = int(start >= 0 and start + width <= len(signal))
            if count % 10_000 == 0:
                logging.info("%s: extracted %d reads", barcode, count)

    metadata["cross_boundary_signal_found"] = found
    metadata["cross_boundary_window_complete"] = complete_window
    metadata["cross_boundary_pre_samples"] = pre_boundary
    metadata["cross_boundary_post_samples"] = post_boundary
    metadata["cross_boundary_index"] = pre_boundary
    np.save(outdir / f"{barcode}.junction_pa_centered.npy", matrix)
    metadata.to_parquet(outdir / f"{barcode}.metadata.parquet", index=False)
    logging.info(
        "%s complete: found=%d/%d complete=%d",
        barcode,
        int(found.sum()),
        len(metadata),
        int(complete_window.sum()),
    )


def main() -> None:
    args = parse_args()
    if args.pre_boundary <= 0 or args.post_boundary <= 0:
        raise ValueError("pre-boundary and post-boundary must be positive")
    args.project_root = args.project_root.resolve()
    args.base_features = args.base_features.resolve()
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)
    os.chdir(args.outdir)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    pod5_dir = (
        args.project_root
        / "MHV_Multiple_sample"
        / "20260306_5bar"
        / "20260306_1823_MN25294_FAZ51928_089bb3a7"
        / "pod5_by_barcode"
    )
    for barcode in BARCODES:
        metadata = pd.read_parquet(
            args.base_features / f"{barcode}.metadata.parquet"
        )
        extract_barcode(
            barcode,
            metadata,
            pod5_dir / f"{barcode}.pod5",
            args.outdir,
            args.pre_boundary,
            args.post_boundary,
        )

    parameters = {
        "pre_boundary_samples": args.pre_boundary,
        "post_boundary_samples": args.post_boundary,
        "boundary_index": args.pre_boundary,
        "signal_units": "pA minus per-read internal-polyA median",
        "sequence_target": "5prime-...A^L U^n-3prime; raw signal order U then A",
    }
    (args.outdir / "cross_boundary_parameters.json").write_text(
        json.dumps(parameters, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
