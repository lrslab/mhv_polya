#!/usr/bin/env python3
"""Validate barcode-figure counts, source tables and reference parameters."""
import argparse
import json
from pathlib import Path

import pandas as pd
from PIL import Image

from reproduce import ROOT


def check_table(frame, rows, columns):
    if len(frame) != rows or list(frame.columns) != columns:
        raise ValueError("Unexpected read-table dimensions or columns")
    if frame["read_id"].isna().any() or frame["read_id"].duplicated().any():
        raise ValueError("Read IDs must be present and unique")


def compare_reads(left, right):
    pd.testing.assert_frame_equal(
        left.sort_values("read_id").reset_index(drop=True),
        right.sort_values("read_id").reset_index(drop=True),
        check_dtype=False, check_exact=True,
    )


def check_calls(frame, cutoff):
    eligible = frame["eligible_terminal_nona_direct"].eq(1)
    fields = ["terminal_nona_signed_difference_pa", "terminal_nona_abs_difference_pa", "terminal_junction_nonA_like"]
    if frame.loc[~eligible, fields].notna().any().any():
        raise ValueError("Uncallable reads contain derived scores or calls")
    if frame.loc[eligible, fields].isna().any().any():
        raise ValueError("Callable reads have missing scores or calls")
    signed = frame.loc[eligible, "terminal_nona_signed_difference_pa"]
    absolute = frame.loc[eligible, "terminal_nona_abs_difference_pa"]
    if not signed.abs().eq(absolute).all():
        raise ValueError("Signed and absolute scores disagree")
    if not absolute.ge(cutoff).eq(frame.loc[eligible, "terminal_junction_nonA_like"]).all():
        raise ValueError("Binary calls disagree with the fixed cutoff")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--reference-results", type=Path)
    p.add_argument("--report-name", default="figure_validation.json")
    a = p.parse_args()
    config = json.loads((ROOT / "config/analysis_config.json").read_text())
    schema = json.loads((ROOT / "schemas/figure_io_schema.json").read_text())
    report = {"status": "RUNNING", "target": "barcode Figure 6A", "checks": []}
    scores = {}
    try:
        for item in schema["tables"]:
            if "internal_path" not in item:
                continue
            columns = [col["name"] for col in item["columns"]]
            frame = pd.read_parquet(a.run / item["internal_path"])
            check_table(frame, item["rows"], columns)
            check_calls(frame, config["fixed_cutoff_pa"])
            if a.reference_results is not None:
                expected = pd.read_parquet(a.reference_results / item["reference_path"], columns=columns)
                compare_reads(frame, expected)
            scores[Path(item["internal_path"]).stem] = frame
            report["checks"].append({"file": item["internal_path"], "rows": len(frame), "columns": len(frame.columns),
                                      "status": "PASS", "external_reference_compared": a.reference_results is not None})
        for item in schema["tables"]:
            if not item["path"].endswith(".tsv.gz"):
                continue
            source = scores.get(Path(item["path"]).name.removesuffix(".tsv.gz"))
            dtypes = {"terminal_junction_nonA_like": "boolean"} if source is None else source.dtypes.to_dict()
            options = {"dtype": dtypes}
            exported = pd.read_csv(a.run / item["path"], sep="\t", float_precision="round_trip", **options)
            check_table(exported, item["rows"], [col["name"] for col in item["columns"]])
            if source is not None:
                compare_reads(exported, source)
            else:
                source = scores["mhv_read_scores"]
                source = source.loc[source["mapping_barcode"].isin(["bc05", "bc06", "bc07", "bc11"])]
                fields = ["read_id", "mapping_barcode", "eligible_terminal_nona_direct", "terminal_junction_nonA_like"]
                compare_reads(exported[fields], source[fields])
            report["checks"].append({"file": item["path"], "rows": len(exported), "status": "PASS"})
        summary = "figure6A/data/figure6A_overall_virus_barcode_summary.tsv"
        pd.testing.assert_frame_equal(
            pd.read_csv(a.run / summary, sep="\t"),
            pd.read_csv(ROOT / "expected/figure6A/figure6A_overall_virus_barcode_summary.tsv", sep="\t"),
            check_dtype=False, check_exact=False, rtol=0, atol=1e-12,
        )
        report["checks"].append({"file": summary, "status": "PASS"})
        params = json.loads((a.run / "results/analysis_parameters.json").read_text())
        host = scores["host_read_scores"]
        reference = host.loc[host["long_host_reference_member"].eq(1), "terminal_minus_stableA_pa"]
        if len(reference) != config["expected_reference_n"] or float(reference.median()) != config["expected_reference_median_pa"]:
            raise ValueError("Pooled host reference differs from the study result")
        for key, expected in [("reference_n", config["expected_reference_n"]),
                              ("reference_median_pa", config["expected_reference_median_pa"]),
                              ("cutoff_pa", config["fixed_cutoff_pa"])]:
            if params[key] != expected:
                raise ValueError(f"Unexpected {key}: {params[key]}")
        report["checks"].append({"file": "results/analysis_parameters.json", "status": "PASS",
                                  "reference_n": params["reference_n"], "reference_median_pa": params["reference_median_pa"]})
        for extension in ("png", "pdf", "svg"):
            rel = f"figure6A/output/figure6A_barcode_level_mut_vs_wt_virus.{extension}"
            path = a.run / rel
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"Missing figure: {path}")
            record = {"file": rel, "bytes": path.stat().st_size, "status": "PASS"}
            if extension == "png":
                with Image.open(path) as picture:
                    record.update({"pixels": list(picture.size), "dpi": list(picture.info.get("dpi", []))})
                    picture.verify()
            report["checks"].append(record)
        report["status"] = "PASS"
    except BaseException as error:
        report.update({"status": "FAILED", "error": f"{type(error).__name__}: {error}"})
        raise
    finally:
        (a.run / a.report_name).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "n_checks": len(report["checks"])}))


if __name__ == "__main__":
    main()
