#!/usr/bin/env python3
"""Call terminal-junction Non-A-like signal using a fixed pA cutoff.

The method compares a terminal dwell with a downstream stable-A window in
the same read, centres the contrast on a pooled long-host reference median,
and applies an absolute cutoff to combine both current directions."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


BARCODES = ("bc04", "bc05", "bc06", "bc07", "bc11")
GROUP_MAP = {
    "bc04": "mock",
    "bc05": "MHV-wt",
    "bc06": "MHV-wt",
    "bc07": "MHV-mut",
    "bc11": "MHV-mut",
}
EXPECTED_VIRUS_COUNTS = {
    "bc04": 1506,
    "bc05": 11405,
    "bc06": 12182,
    "bc07": 17284,
    "bc11": 44029,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--pod5-dir", type=Path)
    parser.add_argument("--sgrna-master", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--anchor-mode",
        choices=("dorado_direct", "warp_fixed", "warp_raw"),
        default="dorado_direct",
    )
    parser.add_argument(
        "--boundary-method",
        choices=("all", "cnn"),
        default="all",
        help="Optionally restrict the fixed Warp boundary to the primary CNN method.",
    )
    parser.add_argument("--cutoff-pa", type=float, default=5.0)
    parser.add_argument("--right-phase-nt", type=float, default=0.5)
    parser.add_argument("--fixed-rna-dwell-samples", type=float, default=31.0)
    parser.add_argument("--warp-dorado-tolerance", type=float, default=20.0)
    parser.add_argument("--min-dorado-polya-nt", type=float, default=10.0)
    parser.add_argument("--min-rna-dwell", type=float, default=20.0)
    parser.add_argument("--max-rna-dwell", type=float, default=50.0)
    parser.add_argument("--min-anchor-span-nt", type=float, default=1.0)
    parser.add_argument("--max-anchor-span-nt", type=float, default=8.0)
    parser.add_argument("--max-baseline-sigma-pa", type=float, default=3.5)
    parser.add_argument("--max-polya-mad-pa", type=float, default=3.0)
    parser.add_argument("--reference-min-polya-nt", type=float, default=100.0)
    parser.add_argument("--reference-min-qscore", type=float, default=12.0)
    return parser.parse_args()


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    p = successes / total
    denominator = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def extract_one_state(
    signal: np.ndarray,
    right_in_array: np.ndarray,
    dwell: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return terminal mean, stable-A mean, contrast, and array-complete flag.

    Terminal window: [right-dwell, right).
    Stable-A window: [right+dwell, right+2*dwell).
    The intervening first A dwell is skipped to reduce boundary-transition
    contamination.  All coordinates use NumPy's round-to-nearest convention.
    """

    n = len(signal)
    terminal = np.full(n, np.nan, dtype=np.float32)
    stable = np.full(n, np.nan, dtype=np.float32)
    contrast = np.full(n, np.nan, dtype=np.float32)
    complete = np.zeros(n, dtype=bool)
    for row in range(n):
        if not np.isfinite(right_in_array[row]) or not np.isfinite(dwell[row]):
            continue
        right = int(round(float(right_in_array[row])))
        width = int(round(float(dwell[row])))
        if width < 1:
            continue
        terminal_lo = right - width
        terminal_hi = right
        stable_lo = right + width
        stable_hi = right + 2 * width
        if terminal_lo < 0 or stable_hi > signal.shape[1]:
            continue
        terminal_segment = np.asarray(signal[row, terminal_lo:terminal_hi], dtype=np.float32)
        stable_segment = np.asarray(signal[row, stable_lo:stable_hi], dtype=np.float32)
        if not np.isfinite(terminal_segment).all() or not np.isfinite(stable_segment).all():
            continue
        terminal[row] = float(terminal_segment.mean())
        stable[row] = float(stable_segment.mean())
        contrast[row] = terminal[row] - stable[row]
        complete[row] = True
    return terminal, stable, contrast, complete


def prepare_barcode_raw(
    barcode: str,
    scored: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Extract the two fixed windows directly from POD5 signal in pA."""

    if args.pod5_dir is None:
        raise ValueError("--pod5-dir is required for --anchor-mode warp_raw")
    import pod5

    pod5_path = args.pod5_dir / f"{barcode}.pod5"
    if not pod5_path.exists():
        raise FileNotFoundError(pod5_path)
    n = len(scored)
    lookup = {str(read_id): index for index, read_id in enumerate(scored["read_id"])}
    terminal = np.full(n, np.nan, dtype=np.float32)
    stable = np.full(n, np.nan, dtype=np.float32)
    contrast = np.full(n, np.nan, dtype=np.float32)
    stable_mad = np.full(n, np.nan, dtype=np.float32)
    signal_found = np.zeros(n, dtype=bool)
    state_complete = np.zeros(n, dtype=bool)
    warp_start_values = scored["polya_start"].to_numpy(np.int64)
    dwell = float(args.fixed_rna_dwell_samples)
    width = int(round(dwell))
    right_offset = int(round(args.right_phase_nt * dwell))
    required_post_boundary = right_offset + 2 * width

    with pod5.Reader(pod5_path) as reader:
        for record in reader.reads(selection=list(lookup)):
            index = lookup.get(str(record.read_id))
            if index is None:
                continue
            signal_found[index] = True
            signal = np.asarray(record.signal_pa, dtype=np.float32)
            warp_start = int(warp_start_values[index])
            right = warp_start + right_offset
            terminal_lo = right - width
            stable_lo = right + width
            stable_hi = right + 2 * width
            if terminal_lo < 0 or stable_hi > len(signal):
                continue
            terminal_segment = signal[terminal_lo:right]
            stable_segment = signal[stable_lo:stable_hi]
            if not np.isfinite(terminal_segment).all() or not np.isfinite(stable_segment).all():
                continue
            terminal[index] = float(terminal_segment.mean())
            stable[index] = float(stable_segment.mean())
            contrast[index] = terminal[index] - stable[index]
            stable_center = float(np.median(stable_segment))
            stable_mad[index] = float(np.median(np.abs(stable_segment - stable_center)))
            state_complete[index] = True

    warp_start = scored["polya_start"].to_numpy(float)
    warp_end = scored["polya_end"].to_numpy(float)
    tail_window_complete = (
        np.isfinite(warp_start)
        & np.isfinite(warp_end)
        & ((warp_end - warp_start) >= required_post_boundary)
    )
    barcode_ok = scored["barcode_mapping_consistent"].fillna(0).eq(1).to_numpy()
    group_ok = scored["mapping_group"].isin(["host", "virus"]).to_numpy()
    boundary_ok = (
        np.ones(n, dtype=bool)
        if args.boundary_method == "all"
        else scored["boundary_method"].eq(args.boundary_method).to_numpy()
    )
    eligible = (
        signal_found
        & group_ok
        & barcode_ok
        & boundary_ok
        & tail_window_complete
        & state_complete
    )

    qc = np.full(n, "eligible", dtype=object)
    qc[~signal_found] = "pod5_signal_not_found"
    prefix = signal_found & group_ok
    qc[prefix & ~barcode_ok] = "barcode_mapping_mismatch"
    prefix &= barcode_ok
    qc[prefix & ~boundary_ok] = "nonprimary_boundary_method"
    prefix &= boundary_ok
    qc[prefix & ~tail_window_complete] = "warp_polya_too_short_for_stableA_window"
    prefix &= tail_window_complete
    qc[prefix & ~state_complete] = "terminal_or_stableA_window_incomplete"

    detected_right = warp_start
    model_right = warp_start + right_offset
    scored["sample_group"] = scored["mapping_barcode"].map(GROUP_MAP)
    scored["terminal_nona_anchor_mode"] = args.anchor_mode
    scored["terminal_nona_coordinate_system"] = "raw_signal_sample"
    scored["terminal_nona_right_phase_nt"] = args.right_phase_nt
    scored["terminal_nona_dorado_warp_delta_samples"] = (
        scored["polya_start"].astype(float) - scored["dorado_polya_start"].astype(float)
    )
    scored["terminal_nona_rna_dwell_samples_per_nt"] = dwell
    scored["terminal_nona_detected_right_array_index"] = detected_right
    scored["terminal_nona_model_right_array_index"] = model_right
    scored["terminal_nona_pa0_array_index"] = scored["dorado_pa_anchor"].to_numpy(float)
    scored["terminal_nona_terminal_left_array_index"] = model_right - width
    scored["terminal_window_mean_centered_pa"] = terminal
    scored["stableA_window_mean_centered_pa"] = stable
    scored["stableA_window_mad_pa"] = stable_mad
    scored["terminal_minus_stableA_pa"] = contrast
    scored["eligible_terminal_nona_direct"] = eligible.astype(np.int8)
    scored["terminal_nona_direct_qc_status"] = qc
    scored["raw_pod5_signal_found"] = signal_found.astype(np.int8)
    scored["warp_tail_window_complete"] = tail_window_complete.astype(np.int8)
    return scored


def prepare_barcode(barcode: str, args: argparse.Namespace) -> pd.DataFrame:
    metadata_path = args.features / f"{barcode}.metadata.parquet"
    signal_path = args.features / f"{barcode}.junction_pa_centered.npy"
    anchor_path = args.anchors / f"{barcode}.dorado_polya_anchors.parquet"
    metadata = pd.read_parquet(metadata_path)
    anchors = pd.read_parquet(anchor_path)
    if metadata["read_id"].duplicated().any() or anchors["read_id"].duplicated().any():
        raise RuntimeError(f"Duplicate read IDs in {barcode} inputs")
    anchors = anchors.set_index("read_id").reindex(metadata["read_id"]).reset_index()
    if args.anchor_mode == "dorado_direct" and anchors["dorado_pa_anchor"].isna().any():
        raise RuntimeError(f"Missing Dorado pa0 anchors for {barcode}")

    keep_anchor = [
        "dorado_polya_estimate_nt",
        "dorado_pa_anchor",
        "dorado_polya_start",
        "dorado_polya_end",
        "dorado_mean_qscore",
    ]
    scored = pd.concat(
        [metadata.reset_index(drop=True), anchors[keep_anchor].reset_index(drop=True)],
        axis=1,
    )
    if args.anchor_mode == "warp_raw":
        return prepare_barcode_raw(barcode, scored, args)
    signal = np.load(signal_path, mmap_mode="r")
    if signal.shape != (len(metadata), 600):
        raise RuntimeError(f"Unexpected signal shape for {barcode}: {signal.shape}")

    valid_tail = (
        scored["dorado_polya_estimate_nt"].ge(args.min_dorado_polya_nt)
        & scored["dorado_polya_start"].ge(0)
        & scored["dorado_polya_end"].gt(scored["dorado_polya_start"])
    )
    delta = scored["polya_start"].astype(float) - scored["dorado_polya_start"].astype(float)
    concordant = valid_tail & delta.abs().le(args.warp_dorado_tolerance)
    with np.errstate(divide="ignore", invalid="ignore"):
        direct_dwell = (
            scored["dorado_polya_end"].to_numpy(float)
            - scored["dorado_polya_start"].to_numpy(float)
        ) / scored["dorado_polya_estimate_nt"].to_numpy(float)

    direct_detected_right = (
        scored["cross_boundary_index"].to_numpy(float)
        + scored["dorado_polya_start"].to_numpy(float)
        - scored["polya_start"].to_numpy(float)
    )
    pa0_in_array = (
        scored["cross_boundary_index"].to_numpy(float)
        + scored["dorado_pa_anchor"].to_numpy(float)
        - scored["polya_start"].to_numpy(float)
    )
    if args.anchor_mode == "dorado_direct":
        dwell = direct_dwell
        detected_right = direct_detected_right
    else:
        dwell = np.full(len(scored), args.fixed_rna_dwell_samples, dtype=float)
        # The signal arrays are centered exactly at WarpDemuX polya_start.
        # Using this common coordinate removes barcode-dependent Dorado pa0/pa1
        # availability from the definition of callability.
        detected_right = scored["cross_boundary_index"].to_numpy(float)
    model_right = detected_right + args.right_phase_nt * dwell
    terminal, stable, contrast, state_complete = extract_one_state(signal, model_right, dwell)

    rounded_right = np.rint(model_right)
    rounded_dwell = np.rint(dwell)
    with np.errstate(invalid="ignore"):
        terminal_left = rounded_right - rounded_dwell
    direct_adapter_safe = (
        np.isfinite(terminal_left)
        & np.isfinite(pa0_in_array)
        & (terminal_left >= pa0_in_array)
    )
    span_nt = (
        scored["dorado_polya_start"].to_numpy(float)
        - scored["dorado_pa_anchor"].to_numpy(float)
    ) / dwell
    direct_valid_anchor = (
        concordant.to_numpy()
        & scored["dorado_pa_anchor"].ge(0).to_numpy()
        & scored["dorado_polya_start"].gt(scored["dorado_pa_anchor"]).to_numpy()
    )
    direct_valid_dwell = np.isfinite(direct_dwell) & (direct_dwell >= args.min_rna_dwell) & (direct_dwell <= args.max_rna_dwell)
    direct_valid_span = np.isfinite(span_nt) & (span_nt >= args.min_anchor_span_nt) & (span_nt <= args.max_anchor_span_nt)
    if args.anchor_mode == "dorado_direct":
        valid_anchor = direct_valid_anchor
        valid_dwell = direct_valid_dwell
        valid_span = direct_valid_span
        adapter_safe = direct_adapter_safe
    else:
        valid_anchor = np.ones(len(scored), dtype=bool)
        valid_dwell = np.isfinite(dwell) & (dwell >= args.min_rna_dwell) & (dwell <= args.max_rna_dwell)
        valid_span = np.ones(len(scored), dtype=bool)
        # warp_fixed deliberately measures the shared mixed DNA|RNA junction;
        # coarse, barcode-dependent Dorado pa0 is retained for QC but not used
        # as an exclusion gate.
        adapter_safe = np.ones(len(scored), dtype=bool)
    base = scored["callable_base"].eq(1).to_numpy()
    barcode_ok = scored["barcode_mapping_consistent"].fillna(0).eq(1).to_numpy()
    cross_complete = scored["cross_boundary_window_complete"].eq(1).to_numpy()
    group_ok = scored["mapping_group"].isin(["host", "virus"]).to_numpy()
    signal_qc = (
        scored["baseline_sigma_pa"].le(args.max_baseline_sigma_pa).to_numpy()
        & scored["polya_mad"].le(args.max_polya_mad_pa).to_numpy()
    )
    boundary_ok = (
        np.ones(len(scored), dtype=bool)
        if args.boundary_method == "all"
        else scored["boundary_method"].eq(args.boundary_method).to_numpy()
    )
    eligible = (
        base
        & barcode_ok
        & cross_complete
        & group_ok
        & signal_qc
        & boundary_ok
        & valid_anchor
        & valid_dwell
        & valid_span
        & adapter_safe
        & state_complete
    )

    qc = np.full(len(scored), "eligible", dtype=object)
    qc[~base] = "insufficient_signal_or_polyA"
    qc[base & ~barcode_ok] = "barcode_mapping_mismatch"
    qc[base & barcode_ok & ~cross_complete] = "cross_boundary_window_incomplete"
    prefix = base & barcode_ok & cross_complete
    qc[prefix & ~signal_qc] = "noisy_internal_polyA_baseline"
    prefix &= signal_qc
    qc[prefix & ~boundary_ok] = "nonprimary_boundary_method"
    prefix &= boundary_ok
    if args.anchor_mode == "dorado_direct":
        qc[prefix & ~concordant.to_numpy()] = "dorado_pa1_not_concordant_with_warp"
        prefix &= concordant.to_numpy()
        qc[prefix & ~valid_anchor] = "invalid_adapter_or_pureA_anchor"
        prefix &= valid_anchor
        qc[prefix & ~valid_dwell] = "invalid_direct_rna_dwell"
        prefix &= valid_dwell
        qc[prefix & ~valid_span] = "adapter_to_pureA_span_out_of_range"
        prefix &= valid_span
        qc[prefix & ~adapter_safe] = "terminal_window_precedes_dorado_pa0"
        prefix &= adapter_safe
    else:
        qc[prefix & ~valid_dwell] = "invalid_fixed_rna_dwell"
        prefix &= valid_dwell
    qc[prefix & ~state_complete] = "terminal_or_stableA_window_incomplete"

    scored["sample_group"] = scored["mapping_barcode"].map(GROUP_MAP)
    scored["terminal_nona_anchor_mode"] = args.anchor_mode
    scored["terminal_nona_right_phase_nt"] = args.right_phase_nt
    scored["terminal_nona_dorado_warp_delta_samples"] = delta
    scored["terminal_nona_rna_dwell_samples_per_nt"] = dwell
    scored["terminal_nona_detected_right_array_index"] = detected_right
    scored["terminal_nona_model_right_array_index"] = model_right
    scored["terminal_nona_pa0_array_index"] = pa0_in_array
    scored["terminal_nona_terminal_left_array_index"] = terminal_left
    scored["terminal_window_mean_centered_pa"] = terminal
    scored["stableA_window_mean_centered_pa"] = stable
    scored["terminal_minus_stableA_pa"] = contrast
    scored["eligible_terminal_nona_direct"] = eligible.astype(np.int8)
    scored["terminal_nona_direct_qc_status"] = qc
    return scored


def summarise(frame: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    grouper: str | list[str] = group_columns[0] if len(group_columns) == 1 else group_columns
    for keys, part in frame.groupby(grouper, dropna=False, observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        eligible = part["eligible_terminal_nona_direct"].eq(1)
        positive = eligible & part["terminal_junction_nonA_like"].eq(True)
        n_total = len(part)
        n_callable = int(eligible.sum())
        n_positive = int(positive.sum())
        lo, hi = wilson_interval(n_positive, n_callable)
        row = {column: value for column, value in zip(group_columns, keys, strict=True)}
        row.update(
            {
                "n_total": n_total,
                "n_callable": n_callable,
                "callable_fraction": n_callable / n_total if n_total else np.nan,
                "n_terminal_junction_nonA_like": n_positive,
                "n_A_junction_like": n_callable - n_positive,
                "terminal_junction_nonA_like_fraction_callable": n_positive / n_callable if n_callable else np.nan,
                "wilson95_low_descriptive": lo,
                "wilson95_high_descriptive": hi,
                "all_read_lower_bound": n_positive / n_total if n_total else np.nan,
                "median_abs_difference_pa_callable": part.loc[eligible, "terminal_nona_abs_difference_pa"].median(),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def length_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    callable_frame = frame.loc[frame["eligible_terminal_nona_direct"].eq(1)].copy()
    for keys, part in callable_frame.groupby(
        ["sample_group", "sgrna_primary_assignment", "terminal_nona_direct_class"],
        dropna=False,
        observed=True,
    ):
        row: dict[str, object] = {
            "sample_group": keys[0],
            "sgrna_primary_assignment": keys[1],
            "terminal_nona_direct_class": keys[2],
            "n_reads": len(part),
        }
        for source, prefix in [("polya_length_nt", "polya_nt"), ("read_nt_len", "read_nt")]:
            values = pd.to_numeric(part[source], errors="coerce").dropna()
            row[f"{prefix}_n"] = len(values)
            row[f"{prefix}_median"] = values.median() if len(values) else np.nan
            row[f"{prefix}_q25"] = values.quantile(0.25) if len(values) else np.nan
            row[f"{prefix}_q75"] = values.quantile(0.75) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def score_all_reads(args: argparse.Namespace):
    """Calculate signal scores before adding transcript annotations."""
    if args.cutoff_pa <= 0:
        raise ValueError("--cutoff-pa must be positive")
    if not 0 <= args.right_phase_nt <= 1:
        raise ValueError("--right-phase-nt must be in [0,1]")
    if args.fixed_rna_dwell_samples <= 0:
        raise ValueError("--fixed-rna-dwell-samples must be positive")
    all_reads = pd.concat(
        [prepare_barcode(barcode, args) for barcode in BARCODES],
        ignore_index=True,
    )
    if all_reads["read_id"].duplicated().any():
        raise RuntimeError("Duplicate read IDs across barcode feature sets")

    if args.anchor_mode == "dorado_direct":
        reference_length_ok = all_reads["dorado_polya_estimate_nt"].ge(
            args.reference_min_polya_nt
        )
        reference_length_source = "Dorado pt"
    else:
        reference_length_ok = all_reads["actual_polya_len_samples"].ge(
            args.reference_min_polya_nt * args.fixed_rna_dwell_samples
        )
        reference_length_source = "WarpDemuX signal span / fixed dwell"
    reference_mask = (
        all_reads["mapping_group"].eq("host")
        & all_reads["eligible_terminal_nona_direct"].eq(1)
        & reference_length_ok
        & all_reads["dorado_mean_qscore"].ge(args.reference_min_qscore)
    )
    reference_values = all_reads.loc[reference_mask, "terminal_minus_stableA_pa"].dropna()
    if len(reference_values) < 1000:
        raise RuntimeError(f"Only {len(reference_values)} long-host reference reads")
    pooled_reference = float(reference_values.median())

    eligible = all_reads["eligible_terminal_nona_direct"].eq(1)
    all_reads["pooled_long_host_A_DNA_reference_pa"] = pooled_reference
    all_reads["terminal_nona_fixed_cutoff_pa"] = args.cutoff_pa
    all_reads["terminal_nona_signed_difference_pa"] = np.nan
    all_reads["terminal_nona_abs_difference_pa"] = np.nan
    signed = all_reads.loc[eligible, "terminal_minus_stableA_pa"] - pooled_reference
    all_reads.loc[eligible, "terminal_nona_signed_difference_pa"] = signed
    all_reads.loc[eligible, "terminal_nona_abs_difference_pa"] = signed.abs()
    all_reads["terminal_junction_nonA_like"] = pd.Series(pd.NA, index=all_reads.index, dtype="boolean")
    all_reads.loc[eligible, "terminal_junction_nonA_like"] = (
        all_reads.loc[eligible, "terminal_nona_abs_difference_pa"].ge(args.cutoff_pa)
    ).to_numpy()
    all_reads["terminal_nona_direct_class"] = "uncallable"
    all_reads.loc[eligible & all_reads["terminal_junction_nonA_like"].eq(True), "terminal_nona_direct_class"] = "terminal_junction_nonA_like"
    all_reads.loc[eligible & all_reads["terminal_junction_nonA_like"].eq(False), "terminal_nona_direct_class"] = "A_junction_like"
    all_reads["long_host_reference_member"] = reference_mask.astype(np.int8)

    return all_reads, reference_mask, pooled_reference, reference_length_source


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    all_reads, reference_mask, pooled_reference, reference_length_source = score_all_reads(args)

    host = all_reads.loc[all_reads["mapping_group"].eq("host")].copy()
    virus = all_reads.loc[all_reads["mapping_group"].eq("virus")].copy()
    observed_counts = virus["mapping_barcode"].value_counts().to_dict()
    if observed_counts != EXPECTED_VIRUS_COUNTS:
        raise RuntimeError(f"Unexpected virus counts: {observed_counts}")

    master = pd.read_parquet(args.sgrna_master)
    if len(master) != sum(EXPECTED_VIRUS_COUNTS.values()) or master["read_id"].duplicated().any():
        raise RuntimeError("Unexpected sgRNA master row count or duplicate read IDs")
    master = master.set_index("read_id").reindex(virus["read_id"]).reset_index()
    if master["mapping_barcode"].to_numpy().tolist() != virus["mapping_barcode"].to_numpy().tolist():
        raise RuntimeError("Barcode mismatch when joining sgRNA master")
    downstream_columns = [
        "sgrna_primary_assignment",
        "sgrna_all_read_identifiability",
        "sgrna_low_confidence_gene_guess",
        "sgrna_assignment_method",
        "sgrna_assignment_confidence",
        "sgrna_primary_high_confidence",
        "polya_length_nt",
        "polya_length_source",
    ]
    for column in downstream_columns:
        virus[column] = master[column].to_numpy()
    if "read_nt_len" not in virus or virus["read_nt_len"].isna().any():
        virus["read_nt_len"] = master["read_nt_len"].to_numpy()
    if virus["read_nt_len"].isna().any() or virus["polya_length_nt"].isna().any():
        raise RuntimeError("Missing read or poly(A) lengths after join")

    group_summary = summarise(virus, ["sample_group"])
    barcode_summary = summarise(virus, ["sample_group", "mapping_barcode"])
    boundary_method_summary = summarise(
        virus,
        ["sample_group", "mapping_barcode", "boundary_method"],
    )
    long_host_reference = all_reads.loc[reference_mask].copy()
    long_host_group_summary = summarise(long_host_reference, ["sample_group"])
    long_host_barcode_summary = summarise(
        long_host_reference,
        ["sample_group", "mapping_barcode"],
    )
    group_comparison = group_summary.merge(
        long_host_group_summary[
            [
                "sample_group",
                "n_callable",
                "n_terminal_junction_nonA_like",
                "terminal_junction_nonA_like_fraction_callable",
            ]
        ].rename(
            columns={
                "n_callable": "n_long_host_reference",
                "n_terminal_junction_nonA_like": "n_long_host_reference_exceeding_cutoff",
                "terminal_junction_nonA_like_fraction_callable": "long_host_reference_fraction_exceeding_cutoff",
            }
        ),
        on="sample_group",
        how="left",
    )
    group_comparison["virus_minus_long_host_reference_fraction"] = (
        group_comparison["terminal_junction_nonA_like_fraction_callable"]
        - group_comparison["long_host_reference_fraction_exceeding_cutoff"]
    )
    barcode_comparison = barcode_summary.merge(
        long_host_barcode_summary[
            [
                "sample_group",
                "mapping_barcode",
                "n_callable",
                "n_terminal_junction_nonA_like",
                "terminal_junction_nonA_like_fraction_callable",
            ]
        ].rename(
            columns={
                "n_callable": "n_long_host_reference",
                "n_terminal_junction_nonA_like": "n_long_host_reference_exceeding_cutoff",
                "terminal_junction_nonA_like_fraction_callable": "long_host_reference_fraction_exceeding_cutoff",
            }
        ),
        on=["sample_group", "mapping_barcode"],
        how="left",
    )
    barcode_comparison["virus_minus_long_host_reference_fraction"] = (
        barcode_comparison["terminal_junction_nonA_like_fraction_callable"]
        - barcode_comparison["long_host_reference_fraction_exceeding_cutoff"]
    )
    high_conf = virus.loc[virus["sgrna_primary_high_confidence"].fillna(0).eq(1)].copy()
    sgrna_group_summary = summarise(high_conf, ["sample_group", "sgrna_primary_assignment"])
    sgrna_barcode_summary = summarise(
        high_conf,
        ["sample_group", "mapping_barcode", "sgrna_primary_assignment"],
    )

    threshold_rows: list[pd.DataFrame] = []
    for threshold in (3.0, 5.0, 7.0):
        temp = virus.copy()
        temp["terminal_junction_nonA_like"] = pd.Series(pd.NA, index=temp.index, dtype="boolean")
        temp.loc[temp["eligible_terminal_nona_direct"].eq(1), "terminal_junction_nonA_like"] = (
            temp.loc[temp["eligible_terminal_nona_direct"].eq(1), "terminal_nona_abs_difference_pa"].ge(threshold)
        ).to_numpy()
        summary = summarise(temp, ["sample_group", "mapping_barcode"])
        summary.insert(0, "cutoff_pa", threshold)
        threshold_rows.append(summary)
    threshold_sensitivity = pd.concat(threshold_rows, ignore_index=True)

    reference_rows: list[dict[str, object]] = []
    for label, part in list(all_reads.loc[reference_mask].groupby("mapping_barcode")) + [("ALL", all_reads.loc[reference_mask])]:
        values = part["terminal_minus_stableA_pa"].dropna()
        center = float(values.median())
        reference_rows.append(
            {
                "mapping_barcode": label,
                "n_long_host_reference": len(values),
                "median_terminal_minus_stableA_pa": center,
                "mad_terminal_minus_stableA_pa": float((values - center).abs().median()),
                "fraction_exceeding_pooled_5pA_cutoff": float((values.sub(pooled_reference).abs() >= args.cutoff_pa).mean()),
            }
        )
    reference_summary = pd.DataFrame(reference_rows)
    lengths = length_summary(high_conf)
    positive_reads = virus.loc[virus["terminal_junction_nonA_like"].eq(True)].copy()

    per_read_parquet = args.outdir / "mhv_reads_terminal_nona_direct.parquet"
    per_read_tsv = args.outdir / "mhv_reads_terminal_nona_direct.tsv.gz"
    positive_tsv = args.outdir / "mhv_terminal_junction_nona_like_reads.tsv.gz"
    virus.to_parquet(per_read_parquet, index=False)
    virus.to_csv(per_read_tsv, sep="\t", index=False, compression="gzip")
    positive_reads.to_csv(positive_tsv, sep="\t", index=False, compression="gzip")
    host.to_parquet(args.outdir / "host_terminal_nona_direct_scores.parquet", index=False)
    write_tsv(group_summary, args.outdir / "mhv_group_nona_direct_summary.tsv")
    write_tsv(barcode_summary, args.outdir / "mhv_barcode_nona_direct_summary.tsv")
    write_tsv(
        boundary_method_summary,
        args.outdir / "mhv_barcode_boundary_method_nona_direct_summary.tsv",
    )
    write_tsv(
        long_host_group_summary,
        args.outdir / "long_host_group_nona_direct_reference_summary.tsv",
    )
    write_tsv(
        long_host_barcode_summary,
        args.outdir / "long_host_barcode_nona_direct_reference_summary.tsv",
    )
    write_tsv(
        group_comparison,
        args.outdir / "mhv_group_vs_long_host_nona_direct_comparison.tsv",
    )
    write_tsv(
        barcode_comparison,
        args.outdir / "mhv_barcode_vs_long_host_nona_direct_comparison.tsv",
    )
    write_tsv(sgrna_group_summary, args.outdir / "mhv_sgrna_group_nona_direct_summary.tsv")
    write_tsv(sgrna_barcode_summary, args.outdir / "mhv_sgrna_barcode_nona_direct_summary.tsv")
    write_tsv(lengths, args.outdir / "mhv_sgrna_nona_direct_length_summary.tsv")
    write_tsv(threshold_sensitivity, args.outdir / "mhv_barcode_nona_direct_threshold_sensitivity.tsv")
    write_tsv(reference_summary, args.outdir / "host_long_polya_pooled_reference_summary.tsv")

    parameters = {
        "endpoint": "terminal_junction_nonA_like",
        "uses_rna_kmer_model": False,
        "uses_trained_classifier": False,
        "uses_empirical_p_or_fdr": False,
        "anchor_mode": args.anchor_mode,
        "boundary_method": args.boundary_method,
        "right_anchor": (
            "Dorado pa1 only; no Warp fallback"
            if args.anchor_mode == "dorado_direct"
            else "WarpDemuX polya_start for every read"
        ),
        "fixed_rna_dwell_samples": (
            args.fixed_rna_dwell_samples
            if args.anchor_mode in {"warp_fixed", "warp_raw"}
            else None
        ),
        "dorado_pa0_pa1_used_as_exclusion_gate": args.anchor_mode == "dorado_direct",
        "right_phase_nt_fixed": args.right_phase_nt,
        "terminal_window": "[right-dwell,right)",
        "stableA_window": "[right+dwell,right+2*dwell)",
        "score": "abs((terminal_mean-stableA_mean)-pooled_long_host_median)",
        "pooled_long_host_reference_pa": pooled_reference,
        "fixed_cutoff_pa": args.cutoff_pa,
        "counts_both_signal_directions": True,
        "reference_definition": {
            "mapping_group": "host",
            "all_barcodes_pooled": list(BARCODES),
            "minimum_polya_nt": args.reference_min_polya_nt,
            "polya_length_source": reference_length_source,
            "minimum_qscore": args.reference_min_qscore,
            "maximum_baseline_sigma_pa": (
                args.max_baseline_sigma_pa
                if args.anchor_mode != "warp_raw"
                else None
            ),
            "maximum_polya_mad_pa": (
                args.max_polya_mad_pa
                if args.anchor_mode != "warp_raw"
                else None
            ),
            "raw_signal_reference_filtering": (
                "finite local windows, boundary scope, barcode consistency, "
                "Warp tail span, tail length, and Q score only"
                if args.anchor_mode == "warp_raw"
                else None
            ),
            "n_reads": int(reference_mask.sum()),
        },
        "group_map": GROUP_MAP,
        "interpretation": "Operational signal-deviation endpoint; not an absolute biochemical non-A fraction without known controls.",
    }
    (args.outdir / "analysis_parameters.json").write_text(
        json.dumps(parameters, indent=2, sort_keys=True) + "\n"
    )

    validation = {
        "status": "PASS",
        "n_virus_rows": len(virus),
        "n_unique_virus_read_ids": int(virus["read_id"].nunique()),
        "virus_counts_by_barcode": {key: int(value) for key, value in observed_counts.items()},
        "n_host_rows": len(host),
        "n_long_host_reference": int(reference_mask.sum()),
        "pooled_reference_pa": pooled_reference,
        "n_virus_callable": int(virus["eligible_terminal_nona_direct"].sum()),
        "n_virus_terminal_junction_nonA_like": int(virus["terminal_junction_nonA_like"].fillna(False).sum()),
        "n_ineligible_with_binary_call": int(
            virus.loc[virus["eligible_terminal_nona_direct"].eq(0), "terminal_junction_nonA_like"].notna().sum()
        ),
        "n_missing_read_length": int(virus["read_nt_len"].isna().sum()),
        "n_missing_polya_length": int(virus["polya_length_nt"].isna().sum()),
        "contains_p_value_columns": any("p_value" in column.lower() or column.lower().endswith("_p") for column in virus.columns),
    }
    if validation["n_ineligible_with_binary_call"] != 0:
        raise RuntimeError("Ineligible reads leaked binary calls")
    (args.outdir / "validation_report.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )

    if args.anchor_mode == "dorado_direct":
        geometry_description = (
            f"Dorado pa1 plus {args.right_phase_nt:.1f} direct RNA dwell; "
            "the terminal window must start at or after Dorado pa0"
        )
    else:
        geometry_description = (
            f"WarpDemuX polya_start plus {args.right_phase_nt:.1f} of one fixed "
            f"{args.fixed_rna_dwell_samples:.1f}-sample dwell for every read; "
            "Dorado pa0/pa1 are not callability gates"
        )
    readme = f"""# Direct fixed-cutoff terminal junction non-A analysis

This result uses no RNA 9-mer model, trained classifier, empirical p value, or
FDR procedure. For each callable read it calculates one fixed terminal-dwell
mean minus one downstream stable-A-dwell mean. The pooled median from long,
clean host poly(A) reads is `{pooled_reference:.4f} pA`. A read is operationally
called `terminal_junction_nonA_like` when its absolute difference from that
pooled reference is at least `{args.cutoff_pa:.1f} pA`; both signal directions
are counted.

The fixed geometry uses {geometry_description}. Boundary-method scope is
`{args.boundary_method}`. All 86,406 MHV-mapped reads are kept; uncallable reads
have no binary call.

This is a signal-deviation endpoint, not a validated biochemical terminal-U or
terminal-non-A fraction. The pooled-barcode baseline is the requested temporary
assumption; barcode-specific host and virus summaries must be inspected because
WT/Mut labels are confounded with barcode.
"""
    (args.outdir / "README.md").write_text(readme)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "python": platform.python_version(),
        "script": str(Path(__file__).resolve()),
        "inputs": {
            "features": str(args.features.resolve()),
            "anchors": str(args.anchors.resolve()),
            "pod5_dir": (
                str(args.pod5_dir.resolve())
                if args.pod5_dir is not None
                else None
            ),
            "sgrna_master": str(args.sgrna_master.resolve()),
        },
    }
    (args.outdir / "result_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
