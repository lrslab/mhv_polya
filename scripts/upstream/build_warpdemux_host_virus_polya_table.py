#!/usr/bin/env python3
import argparse
import csv
import glob
import gzip
import math
import os
import statistics
import subprocess
from collections import Counter, defaultdict


MAPPED_BARCODES = ["bc04", "bc05", "bc06", "bc07", "bc11"]


def normalize_barcode(value):
    value = str(value).strip()
    if value in ("", "-1", "-1.0"):
        return "unclassified"
    try:
        number = int(float(value))
    except ValueError:
        return value
    return f"bc{number:02d}"


def to_float(value):
    value = str(value).strip()
    if value in ("", "NA", "NaN", "nan", "None"):
        return None
    return float(value)


def fmt(value):
    if value is None or value == "":
        return ""
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.6f}"


def percentile(values, q):
    if not values:
        return None
    xs = sorted(values)
    if len(xs) == 1:
        return xs[0]
    pos = (len(xs) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


def load_predictions(run_dir):
    predictions = {}
    duplicates = 0
    conflicts = 0
    for path in sorted(glob.glob(os.path.join(run_dir, "predictions", "barcode_predictions_*.csv.gz"))):
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                read_id = row.get("#read_id", row.get("read_id", "")).strip()
                rec = {
                    "barcode_raw": row["predicted_barcode"].strip(),
                    "barcode": normalize_barcode(row["predicted_barcode"]),
                    "confidence": row.get("confidence_score", "").strip(),
                    "source_file": os.path.basename(path),
                }
                if read_id in predictions:
                    duplicates += 1
                    if predictions[read_id]["barcode"] != rec["barcode"]:
                        conflicts += 1
                    old_conf = to_float(predictions[read_id]["confidence"]) or -1
                    new_conf = to_float(rec["confidence"]) or -1
                    if new_conf > old_conf:
                        predictions[read_id] = rec
                else:
                    predictions[read_id] = rec
    return predictions, duplicates, conflicts


def load_fastq_lengths(base_dir):
    path = os.path.join(base_dir, "Bam", "fastq.read_lengths.merged.tsv")
    lengths = {}
    source_barcodes = {}
    duplicates = 0
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            read_id = row["read_id"].strip()
            length = to_float(row["read_nt_len"])
            if read_id in lengths:
                duplicates += 1
            else:
                lengths[read_id] = length
                source_barcodes[read_id] = row["sample"].strip()
    return path, lengths, source_barcodes, duplicates


def iter_primary_bam(path):
    cmd = ["samtools", "view", "-F", "2304", path]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            flag = int(fields[1])
            yield {
                "read_id": fields[0],
                "mapped": (flag & 4) == 0,
                "reference": fields[2] if (flag & 4) == 0 else "",
                "mapq": int(fields[4]) if (flag & 4) == 0 else None,
            }
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"samtools failed ({rc}): {path}")


def load_mapping(base_dir):
    root = os.path.join(base_dir, "Bam", "aln_host_then_virus")
    candidates = defaultdict(list)
    membership = defaultdict(set)
    duplicate_primary = Counter()
    seen_primary = set()
    for barcode in MAPPED_BARCODES:
        for group in ("host", "virus"):
            path = os.path.join(root, group, f"{barcode}.{group}.sorted.bam")
            for rec in iter_primary_bam(path):
                rid = rec["read_id"]
                membership[rid].add(barcode)
                key = (barcode, group, rid)
                if key in seen_primary:
                    duplicate_primary[(barcode, group)] += 1
                else:
                    seen_primary.add(key)
                if rec["mapped"]:
                    candidates[rid].append({
                        "group": group,
                        "barcode": barcode,
                        "reference": rec["reference"],
                        "mapq": rec["mapq"],
                    })
    return candidates, membership, duplicate_primary


def choose_mapping(read_id, warp_barcode, candidates, membership):
    items = candidates.get(read_id, [])
    same_barcode = [x for x in items if x["barcode"] == warp_barcode]
    pool = same_barcode or items
    if pool:
        # Sequential mapping intends host first, then virus. Prefer host if an
        # unexpected conflict exists, then the highest MAPQ.
        chosen = max(pool, key=lambda x: (x["group"] == "host", x["mapq"] or -1))
        detail = "mapped_same_barcode" if same_barcode else "mapped_different_barcode"
        return chosen, detail
    if warp_barcode in membership.get(read_id, set()):
        return None, "unmapped_host_and_virus"
    if membership.get(read_id):
        return None, "mapping_record_different_barcode_only"
    if warp_barcode not in MAPPED_BARCODES:
        return None, "mapping_bam_unavailable_for_barcode"
    return None, "no_mapping_record"


def load_existing_estimates(run_dir):
    path = os.path.join(run_dir, "boundaries", "polya_nt_estimated.tsv")
    values = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            values[row["read_id"].strip()] = to_float(row["polya_nt_est"])
    return path, values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--warp-run", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--polya-threshold-nt", type=float, default=4.0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    predictions, pred_dups, pred_conflicts = load_predictions(args.warp_run)
    fastq_path, read_lengths, fastq_barcodes, fastq_dups = load_fastq_lengths(args.base)
    mapping_candidates, mapping_membership, mapping_dups = load_mapping(args.base)
    existing_path, existing_estimates = load_existing_estimates(args.warp_run)

    table_path = os.path.join(args.outdir, "all_reads.warpdemux_host_virus_other.polya.tsv")
    summary_path = os.path.join(args.outdir, "warpdemux_polya_summary_by_barcode_and_mapping.tsv")
    qa_path = os.path.join(args.outdir, "warpdemux_merge_quality_report.tsv")
    failed_path = os.path.join(args.outdir, "warpdemux_failed_reason_summary.tsv")

    fields = [
        "read_id", "warpdemux_barcode", "warpdemux_barcode_raw", "warpdemux_barcode_confidence",
        "mapping_group", "mapping_barcode", "barcode_mapping_consistent", "mapping_status_detail",
        "mapping_reference", "mapping_mapq", "warpdemux_polya_len_samples", "read_nt_len",
        "warpdemux_polya_estimate_nt", "warpdemux_polya_estimate_nt_rounded", "polya_gt_4nt",
        "warpdemux_signal_len", "boundary_source_file", "prediction_source_file", "warpdemux_status"
    ]
    summary_values = defaultdict(list)
    summary_counts = Counter()
    seen_boundaries = set()
    boundary_duplicates = 0
    prediction_missing = 0
    length_missing = 0
    mapping_barcode_mismatch = 0
    estimate_compare_n = 0
    estimate_max_abs_diff = 0.0

    with open(table_path, "w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        boundary_files = sorted(glob.glob(os.path.join(args.warp_run, "boundaries", "detected_boundaries_*.csv")))
        for boundary_file in boundary_files:
            with open(boundary_file, newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    read_id = row["read_id"].strip()
                    if read_id in seen_boundaries:
                        boundary_duplicates += 1
                        continue
                    seen_boundaries.add(read_id)

                    pred = predictions.get(read_id)
                    if pred is None:
                        prediction_missing += 1
                        warp_barcode = "unclassified"
                        barcode_raw = ""
                        confidence = ""
                        pred_source = ""
                    else:
                        warp_barcode = pred["barcode"]
                        barcode_raw = pred["barcode_raw"]
                        confidence = pred["confidence"]
                        pred_source = pred["source_file"]

                    selected, detail = choose_mapping(read_id, warp_barcode, mapping_candidates, mapping_membership)
                    if selected is None:
                        group = "other"
                        map_barcode = ""
                        reference = ""
                        mapq = ""
                        consistent = ""
                    else:
                        group = selected["group"]
                        map_barcode = selected["barcode"]
                        reference = selected["reference"]
                        mapq = selected["mapq"]
                        consistent = "1" if map_barcode == warp_barcode else "0"
                        if consistent == "0":
                            mapping_barcode_mismatch += 1

                    signal_len = to_float(row["signal_len"])
                    polya_samples = to_float(row["polya_len"])
                    read_nt_len = read_lengths.get(read_id)
                    if read_nt_len is None:
                        length_missing += 1
                        estimate_nt = None
                    elif signal_len in (None, 0) or polya_samples is None:
                        estimate_nt = None
                    else:
                        estimate_nt = polya_samples * read_nt_len / signal_len

                    old_estimate = existing_estimates.get(read_id)
                    if estimate_nt is not None and old_estimate is not None:
                        estimate_compare_n += 1
                        estimate_max_abs_diff = max(estimate_max_abs_diff, abs(estimate_nt - old_estimate))

                    summary_counts[(warp_barcode, group, "n_reads")] += 1
                    if estimate_nt is not None:
                        summary_counts[(warp_barcode, group, "n_estimated")] += 1
                        if estimate_nt > args.polya_threshold_nt:
                            summary_counts[(warp_barcode, group, "n_gt_threshold")] += 1
                        summary_values[(warp_barcode, group)].append(estimate_nt)

                    writer.writerow({
                        "read_id": read_id,
                        "warpdemux_barcode": warp_barcode,
                        "warpdemux_barcode_raw": barcode_raw,
                        "warpdemux_barcode_confidence": confidence,
                        "mapping_group": group,
                        "mapping_barcode": map_barcode,
                        "barcode_mapping_consistent": consistent,
                        "mapping_status_detail": detail,
                        "mapping_reference": reference,
                        "mapping_mapq": mapq,
                        "warpdemux_polya_len_samples": fmt(polya_samples),
                        "read_nt_len": fmt(read_nt_len),
                        "warpdemux_polya_estimate_nt": fmt(estimate_nt),
                        "warpdemux_polya_estimate_nt_rounded": "" if estimate_nt is None else str(round(estimate_nt)),
                        "polya_gt_4nt": "" if estimate_nt is None else ("1" if estimate_nt > args.polya_threshold_nt else "0"),
                        "warpdemux_signal_len": fmt(signal_len),
                        "boundary_source_file": os.path.basename(boundary_file),
                        "prediction_source_file": pred_source,
                        "warpdemux_status": "success",
                    })

    summary_fields = [
        "warpdemux_barcode", "mapping_group", "n_reads", "n_with_nt_estimate", "pct_with_nt_estimate",
        "n_polya_gt_4nt", "pct_polya_gt_4nt_among_estimated", "mean_nt", "median_nt", "std_nt",
        "min_nt", "p10_nt", "p25_nt", "p75_nt", "p90_nt", "max_nt"
    ]
    barcode_order = ["bc04", "bc05", "bc06", "bc07", "bc11", "bc12", "unclassified"]
    observed_barcodes = sorted({key[0] for key in summary_counts})
    barcode_order += [bc for bc in observed_barcodes if bc not in barcode_order]
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for barcode in barcode_order:
            for group in ("host", "virus", "other"):
                n = summary_counts[(barcode, group, "n_reads")]
                if n == 0:
                    continue
                ne = summary_counts[(barcode, group, "n_estimated")]
                ngt = summary_counts[(barcode, group, "n_gt_threshold")]
                vals = summary_values[(barcode, group)]
                writer.writerow({
                    "warpdemux_barcode": barcode,
                    "mapping_group": group,
                    "n_reads": n,
                    "n_with_nt_estimate": ne,
                    "pct_with_nt_estimate": fmt(100 * ne / n if n else None),
                    "n_polya_gt_4nt": ngt,
                    "pct_polya_gt_4nt_among_estimated": fmt(100 * ngt / ne if ne else None),
                    "mean_nt": fmt(statistics.fmean(vals) if vals else None),
                    "median_nt": fmt(statistics.median(vals) if vals else None),
                    "std_nt": fmt(statistics.stdev(vals) if len(vals) > 1 else None),
                    "min_nt": fmt(min(vals) if vals else None),
                    "p10_nt": fmt(percentile(vals, 0.10)),
                    "p25_nt": fmt(percentile(vals, 0.25)),
                    "p75_nt": fmt(percentile(vals, 0.75)),
                    "p90_nt": fmt(percentile(vals, 0.90)),
                    "max_nt": fmt(max(vals) if vals else None),
                })

    failed_reasons = Counter()
    failed_ids = set()
    failed_duplicates = 0
    for path in sorted(glob.glob(os.path.join(args.warp_run, "failed_reads", "failed_reads_*.csv.gz"))):
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rid = row["read_id"].strip()
                if rid in failed_ids:
                    failed_duplicates += 1
                failed_ids.add(rid)
                failed_reasons[row.get("fail_reason", "unknown").strip() or "unknown"] += 1
    with open(failed_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["fail_reason", "n_rows"])
        for reason, count in failed_reasons.most_common():
            writer.writerow([reason, count])

    qa_rows = [
        ("warpdemux_success_unique_reads", len(seen_boundaries)),
        ("warpdemux_success_duplicate_read_ids", boundary_duplicates),
        ("warpdemux_prediction_unique_reads", len(predictions)),
        ("warpdemux_prediction_duplicate_read_ids", pred_dups),
        ("warpdemux_prediction_conflicting_duplicates", pred_conflicts),
        ("success_reads_missing_prediction", prediction_missing),
        ("success_reads_missing_read_nt_length", length_missing),
        ("success_reads_with_nt_estimate", len(seen_boundaries) - length_missing),
        ("mapping_barcode_mismatch_for_mapped_reads", mapping_barcode_mismatch),
        ("fastq_length_unique_reads", len(read_lengths)),
        ("fastq_length_duplicate_read_ids", fastq_dups),
        ("existing_estimate_comparison_reads", estimate_compare_n),
        ("existing_estimate_max_abs_diff", fmt(estimate_max_abs_diff)),
        ("warpdemux_failed_unique_reads_excluded", len(failed_ids)),
        ("warpdemux_failed_duplicate_read_ids", failed_duplicates),
        ("failed_success_read_id_overlap", len(failed_ids & seen_boundaries)),
        ("source_boundary_pattern", os.path.join(args.warp_run, "boundaries", "detected_boundaries_*.csv")),
        ("source_prediction_pattern", os.path.join(args.warp_run, "predictions", "barcode_predictions_*.csv.gz")),
        ("source_fastq_lengths", fastq_path),
        ("source_existing_estimate_for_validation", existing_path),
    ]
    for (barcode, group), count in sorted(mapping_dups.items()):
        qa_rows.append((f"duplicate_primary_records_{barcode}_{group}", count))
    with open(qa_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "value"])
        writer.writerows(qa_rows)

    warning_path = os.path.join(args.outdir, "SOURCE_AND_UNITS.txt")
    with open(warning_path, "w") as handle:
        handle.write("Source: WarpDemux detected_boundaries_*.csv and barcode_predictions_*.csv.gz\n")
        handle.write("warpdemux_polya_len_samples is the direct WarpDemux output in signal samples.\n")
        handle.write("warpdemux_polya_estimate_nt = polya_len * read_nt_len / signal_len.\n")
        handle.write("Failed WarpDemux reads are excluded from the main table and summarized separately.\n")

    print(table_path)
    print(summary_path)
    print(qa_path)
    print(failed_path)
    print(warning_path)


if __name__ == "__main__":
    main()
