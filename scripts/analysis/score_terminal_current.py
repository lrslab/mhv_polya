#!/usr/bin/env python3
"""Calculate the within-read signal endpoint used by the barcode figure."""
import argparse
import json
from pathlib import Path

from call_terminal_nona_direct_cutoff import score_all_reads

ROOT = Path(__file__).resolve().parents[2]
SCORE_COLUMNS = [
    "read_id", "mapping_barcode", "sample_group", "mapping_group",
    "warpdemux_barcode", "barcode_mapping_consistent", "boundary_method",
    "polya_start", "polya_end", "actual_polya_len_samples",
    "dorado_mean_qscore", "read_nt_len", "raw_pod5_signal_found",
    "warp_tail_window_complete", "terminal_nona_terminal_left_array_index",
    "terminal_nona_model_right_array_index", "terminal_nona_rna_dwell_samples_per_nt",
    "terminal_window_mean_centered_pa", "stableA_window_mean_centered_pa",
    "stableA_window_mad_pa", "terminal_minus_stableA_pa",
    "pooled_long_host_A_DNA_reference_pa", "terminal_nona_fixed_cutoff_pa",
    "terminal_nona_signed_difference_pa", "terminal_nona_abs_difference_pa",
    "eligible_terminal_nona_direct", "terminal_nona_direct_qc_status",
    "terminal_junction_nonA_like", "terminal_nona_direct_class",
    "long_host_reference_member",
]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", type=Path, required=True)
    p.add_argument("--anchors", type=Path, required=True)
    p.add_argument("--pod5-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    config = json.loads((ROOT / "config/analysis_config.json").read_text())
    args = argparse.Namespace(
        features=a.features, anchors=a.anchors, pod5_dir=a.pod5_dir,
        anchor_mode=config["anchor_mode"], boundary_method=config["boundary_method"],
        fixed_rna_dwell_samples=config["fixed_rna_dwell_samples"],
        right_phase_nt=config["right_phase_nt"], cutoff_pa=config["fixed_cutoff_pa"],
        reference_min_polya_nt=config["reference_minimum_polya_nt"],
        reference_min_qscore=config["reference_minimum_qscore"],
    )
    reads, reference_mask, reference, _ = score_all_reads(args)
    scores = a.output / "work/scores"
    results = a.output / "results"
    scores.mkdir(parents=True)
    results.mkdir(exist_ok=True)
    counts = {}
    for group, name in [("host", "host"), ("virus", "mhv")]:
        frame = reads.loc[reads["mapping_group"].eq(group), SCORE_COLUMNS].copy()
        frame.to_parquet(scores / f"{name}_read_scores.parquet", index=False)
        frame.to_csv(results / f"{name}_read_scores.tsv.gz", sep="\t", index=False, compression="gzip")
        counts[name] = {
            "all_reads": len(frame),
            "callable_reads": int(frame["eligible_terminal_nona_direct"].sum()),
            "nonA_like_reads": int(frame["terminal_junction_nonA_like"].fillna(False).sum()),
        }
    parameters = {
        "method": "Within-read terminal minus stable-poly(A), pooled long-host centering",
        "anchor_mode": args.anchor_mode, "boundary_method": args.boundary_method,
        "fixed_rna_dwell_samples": args.fixed_rna_dwell_samples,
        "terminal_window": "[polya_start-16, polya_start+15)",
        "stableA_window": "[polya_start+46, polya_start+77)",
        "signal_units": "pA", "signal_dtype": "float32",
        "reference_minimum_polya_nt": args.reference_min_polya_nt,
        "reference_minimum_qscore": args.reference_min_qscore,
        "reference_n": int(reference_mask.sum()), "reference_median_pa": reference,
        "cutoff_pa": args.cutoff_pa, "counts": counts,
    }
    (results / "analysis_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    print(json.dumps(parameters, indent=2))


if __name__ == "__main__":
    main()
