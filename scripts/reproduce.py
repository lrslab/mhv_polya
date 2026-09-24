#!/usr/bin/env python3
"""Reproduce the MHV within-read terminal-current analysis."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ("bc04", "bc05", "bc06", "bc07", "bc11")
RUN = Path("MHV_Multiple_sample/20260306_5bar/20260306_1823_MN25294_FAZ51928_089bb3a7")
WARP = Path("MHV_Multiple_sample/Warpdemux_polyA/warpdemux_WDX6_rna004_v1_0_20260318_1039_ade68d6a")
PACKAGE = Path("signalu/publication_terminal_nona_package_v1")


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def build_plan(a):
    return build_advanced_plan(a) if a.advanced else build_figure_plan(a)


def build_figure_plan(a):
    project, out = a.project_root, a.output
    code, upstream = ROOT / "scripts/analysis", ROOT / "scripts/upstream"
    required, steps = [], []
    def step(name, script, *args):
        steps.append({"name": name, "argv": [sys.executable, str(script), *map(str, args)]})
    def need(path):
        required.append(str(path))
    for sample in SAMPLES:
        need(project / RUN / "pod5_by_barcode" / f"{sample}.pod5")
    if a.rebuild_inputs:
        view = out / "work/project_view"
        features, anchors = out / "inputs/base_features", out / "inputs/anchors"
        warp = project / WARP
        for folder, pattern in [(warp / "boundaries", "detected_boundaries_*.csv"),
                                (warp / "predictions", "barcode_predictions_*.csv.gz")]:
            for path in sorted(folder.glob(pattern)) or [folder / pattern]:
                need(path)
        need(warp / "boundaries/polya_nt_estimated.tsv")
        need(project / "MHV_Multiple_sample/Bam/fastq.read_lengths.merged.tsv")
        for sample in SAMPLES:
            need(project / "MHV_Multiple_sample/polyA_bam/polya_5pod5_results" / f"{sample}.polya.bam")
            for group in ("host", "virus"):
                need(project / "MHV_Multiple_sample/Bam/aln_host_then_virus" / group / f"{sample}.{group}.sorted.bam")
        step("01_mapping_wdx_metadata", upstream / "build_warpdemux_host_virus_polya_table.py",
             "--base", project / "MHV_Multiple_sample", "--warp-run", warp,
             "--outdir", view / "polyA_host_virus_by_barcode_warpdemux")
        step("02_signal_metadata", upstream / "extract_terminal_u_features.py",
             "--project-root", view, "--outdir", features)
        step("03_dorado_metadata", upstream / "extract_dorado_polya_anchors.py",
             "--bam-dir", project / "MHV_Multiple_sample/polyA_bam/polya_5pod5_results",
             "--boundary-dir", warp / "boundaries", "--outdir", anchors)
    else:
        features, anchors = a.inputs / "features", a.inputs / "anchors"
        for sample in SAMPLES:
            need(features / f"{sample}.metadata.parquet")
            need(anchors / f"{sample}.dorado_polya_anchors.parquet")
    step("04_terminal_current", code / "score_terminal_current.py", "--features", features,
         "--anchors", anchors, "--pod5-dir", project / RUN / "pod5_by_barcode", "--output", out)
    step("05_barcode_figure", code / "build_figure6A_overall_mut_vs_wt.py",
         "--source-parquet", out / "work/scores/mhv_read_scores.parquet",
         "--output-root", out / "figure6A", "--figure-mode", "barcode")
    reference_args = []
    if a.reference_results is not None:
        for name in ("host_terminal_nona_direct_scores.parquet", "mhv_reads_terminal_nona_direct.parquet"):
            need(a.reference_results / "final_primary" / name)
        reference_args = ["--reference-results", a.reference_results]
    step("06_validate_figure", ROOT / "scripts/validate_figure.py", "--run", out, *reference_args)
    return {"target": "publication_terminal_nona_results_v1", "analysis": "barcode_figure",
            "rebuild_inputs": a.rebuild_inputs, "project_root": str(project), "output": str(out),
            "required_inputs": sorted(set(required)), "steps": steps}


def build_advanced_plan(a):
    project, out = a.project_root, a.output
    code, upstream = ROOT / "scripts/analysis", ROOT / "scripts/upstream"
    input_root = a.inputs or project / PACKAGE / "inputs"
    reference_root = a.reference_results
    steps, required = [], []
    def step(name, script, *args):
        steps.append({"name": name, "argv": [sys.executable, str(script), *map(str, args)]})
    def need(path):
        required.append(str(path))
    for sample in SAMPLES:
        need(project / RUN / "pod5_by_barcode" / f"{sample}.pod5")
    if a.rebuild_inputs:
        view = out / "work/project_view"
        feature_root, anchors = out / "inputs/features", out / "inputs/anchors"
        master = out / "inputs/sgrna_master.parquet"
        warp = project / WARP
        for folder, pattern in [(warp / "boundaries", "detected_boundaries_*.csv"), (warp / "predictions", "barcode_predictions_*.csv.gz")]:
            files = sorted(folder.glob(pattern))
            if files:
                for path in files: need(path)
            else: need(folder / pattern)
        need(warp / "boundaries/polya_nt_estimated.tsv")
        need(project / "MHV_Multiple_sample/Bam/fastq.read_lengths.merged.tsv")
        for sample in SAMPLES:
            need(project / "MHV_Multiple_sample/polyA_bam/polya_5pod5_results" / f"{sample}.polya.bam")
            for origin in ("host", "virus"):
                bam = project / "MHV_Multiple_sample/Bam/aln_host_then_virus" / origin / f"{sample}.{origin}.sorted.bam"
                need(bam)
                if origin == "virus": need(Path(str(bam) + ".bai"))
        step("01_mapping_wdx_metadata", upstream / "build_warpdemux_host_virus_polya_table.py",
            "--base", project / "MHV_Multiple_sample", "--warp-run", warp,
            "--outdir", view / "polyA_host_virus_by_barcode_warpdemux")
        step("02_base_features", upstream / "extract_terminal_u_features.py", "--project-root", view, "--outdir", out / "inputs/base_features")
        step("03_cross_boundary_features", upstream / "extract_cross_boundary_features.py", "--project-root", view,
            "--base-features", out / "inputs/base_features", "--outdir", feature_root, "--pre-boundary", 240, "--post-boundary", 360)
        step("04_dorado_metadata", upstream / "extract_dorado_polya_anchors.py", "--bam-dir", project / "MHV_Multiple_sample/polyA_bam/polya_5pod5_results",
            "--boundary-dir", warp / "boundaries", "--outdir", anchors)
        step("05_transcript_metadata", upstream / "build_sgrna_master.py", "--features", feature_root, "--anchors", anchors,
            "--bam-dir", project / "MHV_Multiple_sample/Bam/aln_host_then_virus/virus", "--output", master)
    else:
        feature_root, anchors = input_root / "features", input_root / "anchors"
        master = input_root / "sgrna/mhv_sgrna_assignment_length_master.parquet"
        for sample in SAMPLES:
            need(feature_root / f"{sample}.metadata.parquet")
            need(anchors / f"{sample}.dorado_polya_anchors.parquet")
        need(master)
    primary, split = out / "results/final_primary", out / "results/final_host_virus_sgrna"
    if reference_root is not None:
        for name in ("mhv_reads_terminal_nona_direct.parquet", "host_terminal_nona_direct_scores.parquet",
                     "mhv_reads_terminal_nona_direct.tsv.gz", "mhv_terminal_junction_nona_like_reads.tsv.gz"):
            need(reference_root / "final_primary" / name)
        need(reference_root / "final_host_virus_sgrna/virus_high_conf_sgrna_nonA_like_reads.tsv.gz")
        for folder in ("final_primary", "final_host_virus_sgrna"):
            for path in sorted((reference_root / folder).glob("*.tsv")): need(path)
    step("06_self_contrast_5pa", code / "call_terminal_nona_direct_cutoff.py", "--features", feature_root,
        "--anchors", anchors, "--pod5-dir", project / RUN / "pod5_by_barcode", "--sgrna-master", master,
        "--outdir", primary, "--anchor-mode", "warp_raw", "--boundary-method", "cnn",
        "--fixed-rna-dwell-samples", "30.7692307692", "--cutoff-pa", 5, "--right-phase-nt", 0.5)
    step("07_host_virus_sgrna", code / "summarize_host_virus_sgrna_direct.py", "--source-result", primary, "--outdir", split)
    reference_args = [] if reference_root is None else ["--reference-caller", reference_root / "final_primary",
                                                       "--reference-split", reference_root / "final_host_virus_sgrna"]
    step("08_original_validator", code / "validate_reproduction.py", "--caller-dir", primary, "--split-dir", split,
        "--output", out / "analysis_validation.json", *reference_args)
    step("09_publication_figures", code / "plot_publication_figures.py", "--package-root", out)
    step("10_figure6A", code / "build_figure6A_overall_mut_vs_wt.py", "--source-parquet", primary / "mhv_reads_terminal_nona_direct.parquet",
        "--output-root", out / "figure6A_overall_mut_vs_wt", "--figure-mode", "all")
    target_args = [] if reference_root is None else ["--reference-results", reference_root]
    step("11_local_target_comparison", ROOT / "scripts/validate_target.py", "--run", out, *target_args)
    return {"target": "publication_terminal_nona_results_v1", "analysis": "advanced", "rebuild_inputs": a.rebuild_inputs,
        "project_root": str(project), "output": str(out), "required_inputs": sorted(set(required)), "steps": steps}


def preflight(plan):
    errors = []
    output = Path(plan["output"])
    if output.exists() or output.is_symlink(): errors.append(f"Output already exists: {output}")
    for name in plan["required_inputs"]:
        path = Path(name)
        if not path.is_file() or path.stat().st_size == 0: errors.append(f"Missing/empty input: {path}")
    for name in ("numpy", "pandas", "pyarrow", "pod5", "pysam", "scipy", "matplotlib"):
        try: importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: errors.append(f"Missing Python package: {name}")
    if plan["rebuild_inputs"]:
        import shutil
        if not shutil.which("samtools"): errors.append("samtools required to read archived mapping BAMs")
    return errors


def execute(plan):
    out = Path(plan["output"])
    out.mkdir(parents=True, exist_ok=False)
    for name in ("logs", "work", "inputs", "results"):
        (out / name).mkdir()
    if plan["rebuild_inputs"]:
        view = out / "work/project_view"
        view.mkdir()
        (view / "MHV_Multiple_sample").symlink_to(Path(plan["project_root"]) / "MHV_Multiple_sample", target_is_directory=True)
    report = {"status": "RUNNING", "started_utc": now(), "plan": plan, "inputs": [], "steps": [], "packages": {}}
    for name in ("numpy", "pandas", "pyarrow", "pod5", "pysam", "scipy", "matplotlib"):
        report["packages"][name] = importlib.metadata.version(name)
    write_json(out / "run_manifest.json", report)
    try:
        for name in plan["required_inputs"]:
            path = Path(name)
            report["inputs"].append({"path": name, "bytes": path.stat().st_size})
        write_json(out / "run_manifest.json", report)
        env = os.environ.copy()
        env.update({"PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0", "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(out / "work/matplotlib")})
        for step in plan["steps"]:
            print(f"[{now()}] {step['name']}: {shlex.join(step['argv'])}", flush=True)
            item = {**step, "started_utc": now()}
            report["steps"].append(item)
            write_json(out / "run_manifest.json", report)
            with (out / "logs" / f"{step['name']}.log").open("xb") as log:
                result = subprocess.run(step["argv"], cwd=out / "work", env=env, stdout=log, stderr=subprocess.STDOUT)
            item.update({"finished_utc": now(), "returncode": result.returncode})
            if result.returncode: raise RuntimeError(f"{step['name']} failed; see {out / 'logs'}")
        report["status"] = "PASS"
    except BaseException as error:
        report.update({"status": "FAILED", "error": f"{type(error).__name__}: {error}"})
        raise
    finally:
        report["finished_utc"] = now()
        write_json(out / "run_manifest.json", report)
    print(f"PASS: {out}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--inputs", type=Path, help="Alternative packaged inputs directory")
    p.add_argument("--reference-results", type=Path, help="Compare read-level values directly with archived result tables")
    p.add_argument("--advanced", action="store_true", help="Include transcript assignments, sgRNA summaries and Figures 1–4")
    p.add_argument("--rebuild-inputs", action="store_true", help="Rebuild intermediates from the archived source inputs (the default)")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--plan", action="store_true")
    mode.add_argument("--execute", action="store_true")
    a = p.parse_args(argv)
    for key in ("project_root", "inputs", "reference_results"):
        if getattr(a, key): setattr(a, key, getattr(a, key).expanduser().resolve())
    a.output = a.output.expanduser().absolute()
    if a.rebuild_inputs and a.inputs: p.error("Choose --rebuild-inputs or --inputs")
    a.rebuild_inputs = a.rebuild_inputs or a.inputs is None
    plan = build_plan(a)
    if a.plan:
        print(json.dumps(plan, indent=2)); return 0
    errors = preflight(plan)
    if errors:
        print("Preflight FAILED:\n" + "\n".join(errors), file=sys.stderr); return 2
    if not a.execute:
        print(f"Preflight PASS: {len(plan['steps'])} stages. Add --execute to run."); return 0
    execute(plan); return 0


if __name__ == "__main__":
    raise SystemExit(main())
