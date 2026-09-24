#!/usr/bin/env python3
"""Compare advanced result summaries and read-table structure with the study."""
import argparse
import importlib.util
import json
from pathlib import Path

import pandas as pd

from reproduce import ROOT
from validate_figure import check_table, check_calls, compare_reads

spec = importlib.util.spec_from_file_location("original_validation", ROOT / "scripts/analysis/validate_reproduction.py")
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--reference-results", type=Path, help="Compare read-level values with archived result tables")
    p.add_argument("--report-name", default="local_target_validation.json")
    a = p.parse_args()
    report = {"target": "publication_terminal_nona_results_v1", "status": "RUNNING", "checks": []}
    try:
        schema = json.loads((ROOT / "schemas/observed_io_schema.json").read_text())
        for item in schema["tables"]:
            rel = Path(item["path"])
            if item["scope"] != "output" or not (rel.suffix == ".parquet" or str(rel).endswith(".tsv.gz")):
                continue
            frame = pd.read_parquet(a.run / rel) if rel.suffix == ".parquet" else pd.read_csv(a.run / rel, sep="\t", float_precision="round_trip")
            check_table(frame, item["rows"], [col["name"] for col in item["columns"]])
            if rel.suffix == ".parquet":
                check_calls(frame, 5.0)
            compared = False
            if a.reference_results is not None:
                ref = (a.reference_results / Path(*rel.parts[1:]) if rel.parts[0] == "results"
                       else a.reference_results.parent / rel)
                if rel.suffix == ".parquet":
                    compare_reads(frame, pd.read_parquet(ref))
                    compared = True
                elif ref.is_file():
                    compare_reads(frame, pd.read_csv(ref, sep="\t", float_precision="round_trip"))
                    compared = True
            report["checks"].append({"file": str(rel), "rows": len(frame), "columns": len(frame.columns),
                                      "status": "PASS", "external_reference_compared": compared})
        locations = [("final_primary", "results/final_primary"), ("final_host_virus_sgrna", "results/final_host_virus_sgrna"),
                     ("figure_source_data", "figures/source_data"), ("figure6A", "figure6A_overall_mut_vs_wt/data")]
        for expected_dir, run_dir in locations:
            for expected in sorted((ROOT / "expected" / expected_dir).glob("*.tsv")):
                observed = a.run / run_dir / expected.name
                if not observed.is_file() or not original.tables_equal(observed, expected):
                    raise RuntimeError(f"Target TSV differs or is absent: {run_dir}/{expected.name}")
                report["checks"].append({"file": f"{run_dir}/{expected.name}", "status": "PASS"})
        report["status"] = "PASS"
    except BaseException as error:
        report.update({"status": "FAILED", "error": f"{type(error).__name__}: {error}"})
        raise
    finally:
        (a.run / a.report_name).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "n_checks": len(report["checks"])}))


if __name__ == "__main__":
    main()
