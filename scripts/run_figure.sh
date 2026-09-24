#!/usr/bin/env bash
# Usage: bash scripts/run_figure.sh /path/to/MHV_polyA /path/to/barcode_figure_run
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: bash scripts/run_figure.sh INPUT_DIRECTORY OUTPUT_DIRECTORY" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
project_root="$(cd "$1" && pwd)"
mkdir -p "$(dirname "$2")"
mkdir "$2"
output="$(cd "$2" && pwd)"
python_bin="${PYTHON:-python}"

base="${project_root}/MHV_Multiple_sample"
warp="${base}/Warpdemux_polyA/warpdemux_WDX6_rna004_v1_0_20260318_1039_ade68d6a"
pod5="${base}/20260306_5bar/20260306_1823_MN25294_FAZ51928_089bb3a7/pod5_by_barcode"
view="${output}/work/project_view"
mkdir -p "$view" "${output}/inputs" "${output}/results" "${output}/logs"
ln -s "$base" "${view}/MHV_Multiple_sample"
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 MPLBACKEND=Agg
export MPLCONFIGDIR="${output}/work/matplotlib"
cd "${output}/work"

# 1. Join mapping assignments, WDX boundaries and read metadata.
echo "[1/6] Mapping and WDX metadata"
"$python_bin" "${repo_root}/scripts/upstream/build_warpdemux_host_virus_polya_table.py" \
  --base "$base" --warp-run "$warp" \
  --outdir "${view}/polyA_host_virus_by_barcode_warpdemux" \
  > "${output}/logs/01_mapping_wdx_metadata.log" 2>&1

# 2. Extract signal metadata and signal-clipped poly(A) spans.
echo "[2/6] Signal metadata"
"$python_bin" "${repo_root}/scripts/upstream/extract_terminal_u_features.py" \
  --project-root "$view" --outdir "${output}/inputs/base_features" \
  > "${output}/logs/02_signal_metadata.log" 2>&1

# 3. Extract Dorado read quality for host-reference selection.
echo "[3/6] Dorado quality and metadata"
"$python_bin" "${repo_root}/scripts/upstream/extract_dorado_polya_anchors.py" \
  --bam-dir "${base}/polyA_bam/polya_5pod5_results" \
  --boundary-dir "${warp}/boundaries" --outdir "${output}/inputs/anchors" \
  > "${output}/logs/03_dorado_metadata.log" 2>&1

# 4. Compare terminal and stable-poly(A) current within each read.
echo "[4/6] Within-read current analysis"
"$python_bin" "${repo_root}/scripts/analysis/score_terminal_current.py" \
  --features "${output}/inputs/base_features" --anchors "${output}/inputs/anchors" \
  --pod5-dir "$pod5" --output "$output" \
  > "${output}/logs/04_terminal_current.log" 2>&1

# 5. Draw the barcode-level figure and export its source tables.
echo "[5/6] Barcode figure"
"$python_bin" "${repo_root}/scripts/analysis/build_figure6A_overall_mut_vs_wt.py" \
  --source-parquet "${output}/work/scores/mhv_read_scores.parquet" \
  --output-root "${output}/figure6A" --figure-mode barcode \
  > "${output}/logs/05_barcode_figure.log" 2>&1

# 6. Compare the outputs with the expected study results.
echo "[6/6] Result validation"
"$python_bin" "${repo_root}/scripts/validate_figure.py" --run "$output" \
  > "${output}/logs/06_validate_figure.log" 2>&1
echo "PASS: ${output}/figure6A/output/"
