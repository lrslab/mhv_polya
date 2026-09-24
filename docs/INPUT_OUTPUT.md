# Input and output specification: barcode Figure 6A

The figure command is `bash scripts/run_figure.sh P O`, where `P` contains the archived study inputs and `O` is a new output directory. [inputs.tsv](../config/inputs.tsv) supplies exact paths and study file sizes. The optional sgRNA workflow has its own [advanced specification](ADVANCED_INPUT_OUTPUT.md).

## Required inputs: 36 files

All paths are relative to `P`. Expand `bcXX` to `bc04`, `bc05`, `bc06`, `bc07`, `bc11`, and `{0..6}` to the seven WDX shard numbers.

```text
S = MHV_Multiple_sample
D = S/20260306_5bar/20260306_1823_MN25294_FAZ51928_089bb3a7/pod5_by_barcode
B = S/Bam/aln_host_then_virus
W = S/Warpdemux_polyA/warpdemux_WDX6_rna004_v1_0_20260318_1039_ade68d6a
```

| Path | Count | Content used |
|---|---:|---|
| `D/bcXX.pod5` | 5 | Original read UUID, calibrated raw current and signal length |
| `B/host/bcXX.host.sorted.bam` | 5 | Read name, primary alignment flag, reference and MAPQ |
| `B/virus/bcXX.virus.sorted.bam` | 5 | Read name, primary alignment flag, reference and MAPQ |
| `S/polyA_bam/polya_5pod5_results/bcXX.polya.bam` | 5 | Dorado quality and diagnostic metadata |
| `W/boundaries/detected_boundaries_{0..6}.csv` | 7 | WDX signal boundaries and detector metadata |
| `W/predictions/barcode_predictions_{0..6}.csv.gz` | 7 | WDX barcode predictions |
| `S/Bam/fastq.read_lengths.merged.tsv` | 1 | Basecalled read lengths |
| `W/boundaries/polya_nt_estimated.tsv` | 1 | WDX length estimates for metadata QA |

The processing boundary is demultiplexed raw signal plus archived preprocessing outputs. [Reproducibility](REPRODUCIBILITY.md) records the earlier basecalling, demultiplexing and alignment steps. Public data accession is pending.

### Identifiers, groups and coordinates

All joins use the original read UUID: POD5 `read_id`, BAM query name and text-table `read_id`. Mapping barcodes identify the source library. Numeric WDX predictions are normalized to `bc04`, etc. Group labels are `mock` (bc04), `MHV-wt` (bc05/bc06), and `MHV-mut` (bc07/bc11).

WDX positions are zero-based indices into the original POD5 signal. Windows are half-open. The sequencing rate is 4,000 samples/s; the fixed nominal dwell is 30.7692307692 samples/nt. Current is pA, signal lengths are samples, and basecalled read lengths are nt.

### POD5 and BAM

`ReadRecord.signal_pa` supplies calibrated current as float32. The mapping join scans BAM records with `samtools view -F 2304`, excluding secondary and supplementary records. Unmapped primary records contribute barcode membership. Among multiple mapped candidates, the historical join prefers the WDX-matching barcode, then host over virus, then larger MAPQ. This defines host/MHV membership independently of transcript assignment.

Dorado BAMs are scanned sequentially; read IDs are unique within and across the five files. The extractor consumes:

| Tag | Use and missing-value behavior |
|---|---|
| `qs` | Mean read quality. A missing value excludes a host read from the pooled reference. |
| `pt` | Poly(A)-length estimate in nt, retained as intermediate metadata. |
| `pa` | Signal-coordinate array: search anchor, primary start/end, optional secondary start/end. Missing/short arrays are padded with nulls. |
| `ts`, `ns` | Optional trimmed start and signal end; nullable metadata. |
| `sm`, `sd` | Optional scaling midpoint and dispersion; nullable metadata. |

The primary endpoint uses WDX coordinates; Dorado `qs` controls host-reference selection.

### WDX boundary columns

The boundary and Dorado extractors consume these names from the 55-column CSV exports:

```text
read_id, signal_len, adapter_start, adapter_end,
polya_start, polya_end, polya_len,
polya_mean, polya_std, polya_med, polya_mad,
llr_adapter_end, llr_polya_end,
cnn_adapter_end, cnn_polya_end,
start_peak_adapter_end, start_peak_polya_end,
adapter_dt_med, adapter_dt_mad
```

`read_id` is a string. Signal lengths and boundary positions are numeric sample counts. Selected `polya_start`/`polya_end` positions are finite integral values. Alternative detector coordinates can be missing. `boundary_method` is `cnn` when both CNN coordinates are present, otherwise `fallback`. WDX summary-current columns are upstream metadata; endpoint window means are recalculated from POD5.

### Prediction and length tables

| Table | Columns consumed | Format |
|---|---|---|
| WDX predictions | `#read_id` or `read_id`, `predicted_barcode`; optional `confidence_score` | Gzip-compressed CSV |
| FASTQ lengths | `sample`, `read_id`, `read_nt_len` | TSV; 226,920 data rows |
| WDX length estimates | `read_id`, `polya_nt_est` | TSV; 268,378 data rows |

Prediction duplicates retain the larger confidence; read-length duplicates retain the first record. Boundary IDs are unique across shards. Optional WDX failed-read exports contribute QA counts under `work/`; the study reconstruction produced identical analysis metadata using the declared inputs alone.

## Primary outputs

Paths are relative to `O`. TSVs use UTF-8, tabs and a header row. Empty cells are missing values; booleans use `True`/`False`. Fractions use a 0–1 scale, and figure labels convert them to percent.

| Path | Rows × columns / format | Population or purpose |
|---|---|---|
| `results/host_read_scores.tsv.gz` | 82,709 × 30 | All selected host reads, including reference membership |
| `results/mhv_read_scores.tsv.gz` | 86,406 × 30 | All selected MHV reads, including mock and uncallable reads |
| `results/analysis_parameters.json` | JSON | Geometry, threshold, reference and population counts |
| `figure6A/data/figure6A_overall_virus_read_level.tsv.gz` | 84,900 × 19 | MHV reads from bc05/bc06/bc07/bc11, including uncallable reads |
| `figure6A/data/figure6A_overall_virus_barcode_summary.tsv` | 4 × 13 | One row per infected barcode |
| `figure6A/output/figure6A_barcode_level_mut_vs_wt_virus.pdf` | Vector PDF | Publication figure |
| Same stem with `.svg` | Vector SVG | Editable vector figure |
| Same stem with `.png` | 400 dpi PNG | Raster figure |
| `figure_validation.json` | JSON | Ten numerical and output checks |

The [primary schema](../schemas/figure_io_schema.json) records ordered columns, types and missing-value counts from the validated run. Internal score Parquets mirror the 30-column exports and retain the numeric values for direct comparison.

### All 30 score fields

The following order applies to both host and MHV exports.

| Field | Type / unit | Meaning |
|---|---|---|
| `read_id` | String | Original read UUID; unique within each table |
| `mapping_barcode` | Category | Source barcode |
| `sample_group` | Category | `mock`, `MHV-wt`, `MHV-mut` |
| `mapping_group` | Category | `host` or `virus` |
| `warpdemux_barcode` | Category | Normalized WDX assignment |
| `barcode_mapping_consistent` | Numeric 0/1 | WDX/mapping barcode agreement |
| `boundary_method` | Category | `cnn` or `fallback` |
| `polya_start` | Samples | W, WDX poly(A) start |
| `polya_end` | Samples | WDX poly(A) end |
| `actual_polya_len_samples` | Samples | Signal-clipped span: `max(0, min(polya_end, signal_length) - polya_start)` |
| `dorado_mean_qscore` | Numeric, nullable | Dorado `qs` |
| `read_nt_len` | nt | Basecalled read length |
| `raw_pod5_signal_found` | Integer 0/1 | Signal was available |
| `warp_tail_window_complete` | Integer 0/1 | Detected span `polya_end - polya_start >= 77` |
| `terminal_nona_terminal_left_array_index` | Samples | W − 16, inclusive terminal-window left edge |
| `terminal_nona_model_right_array_index` | Samples | W + 15, exclusive terminal-window right edge |
| `terminal_nona_rna_dwell_samples_per_nt` | Samples/nt | Fixed 30.7692307692 |
| `terminal_window_mean_centered_pa` | pA | Absolute calibrated terminal mean T |
| `stableA_window_mean_centered_pa` | pA | Absolute calibrated stable-poly(A) mean A |
| `stableA_window_mad_pa` | pA | Median absolute deviation within the stable-A window |
| `terminal_minus_stableA_pa` | pA | Within-read contrast x = T − A |
| `pooled_long_host_A_DNA_reference_pa` | pA | Pooled host reference x_A, repeated on every read |
| `terminal_nona_fixed_cutoff_pa` | pA | 5.0 |
| `terminal_nona_signed_difference_pa` | pA, nullable | x − x_A for callable reads |
| `terminal_nona_abs_difference_pa` | pA, nullable | Absolute deviation D for callable reads |
| `eligible_terminal_nona_direct` | Integer 0/1 | Combined callability indicator |
| `terminal_nona_direct_qc_status` | Category | Callability or first applicable QC reason |
| `terminal_junction_nonA_like` | Nullable Boolean | True if callable and D ≥5; False if callable and D <5; missing if uncallable |
| `terminal_nona_direct_class` | Category | `terminal_junction_nonA_like`, `A_junction_like`, `uncallable` |
| `long_host_reference_member` | Integer 0/1 | Membership in the 26,997-read host reference; always 0 for MHV |

The historical names ending in `mean_centered_pa` contain **absolute calibrated means**. Centering occurs when the within-read contrast is compared with the pooled host reference.

QC labels are `eligible`, `pod5_signal_not_found`, `barcode_mapping_mismatch`, `nonprimary_boundary_method`, `warp_polya_too_short_for_stableA_window`, and `terminal_or_stableA_window_incomplete`, in that order of precedence. Raw means may be finite on a read that fails another eligibility criterion; its deviation score and call remain missing.

### Figure source fields and denominators

The 19-column figure read table contains `read_id`, `mapping_barcode`, `sample_group`, `display_group`, `mapping_group`, `barcode_mapping_consistent`, `boundary_method`, `eligible_terminal_nona_direct`, `terminal_nona_direct_qc_status`, the two window means, the within-read contrast, the pooled reference, the cutoff, the signed and absolute deviations, the binary call, `terminal_nona_direct_class`, and `figure_signal_class`. `display_group` is WT or Mut; `figure_signal_class` mirrors the callability/classification for plotting. Floating-point values use eight decimal places.

The barcode summary contains three grouping columns (`sample_group`, `display_group`, `mapping_barcode`), four counts (`n_total_mhv`, `n_callable`, `n_uncallable`, `n_terminal_junction_nonA_like`), `n_A_junction_like`, `callable_fraction`, the two class fractions, and the lower/upper Wilson 95% limits for the Non-A-like fraction. Its ordered header is retained in the expected TSV and primary schema.

For each barcode, Non-A-like fraction = positive callable reads / callable reads; A-junction-like fraction = 1 − Non-A-like fraction. Uncallable reads contribute to `n_total_mhv` and `n_uncallable`. Wilson limits describe sampling among reads. The displayed figure is a descriptive barcode comparison.

## Intermediates and execution options

`inputs/base_features/` contains five 39-column metadata Parquets and ten 360-sample NPY arrays produced by the original extractor. `inputs/anchors/` contains five 22-column Dorado tables and a summary TSV. The caller calculates both endpoint windows directly from POD5. `work/` contains the metadata join, two 30-column score Parquets, a read-only project-view symlink and plotting cache. `logs/` contains six stage logs.

The optional `scripts/run_pipeline.sh --project-root P --output O` runner adds a `run_manifest.json` with input paths, commands and package versions. It accepts these additional options: `--plan` prints commands and required paths. Omitting `--execute` runs preflight. `--rebuild-inputs` explicitly selects the default source-input workflow. `--inputs C` reuses five `C/features/bcXX.metadata.parquet` and five `C/anchors/bcXX.dorado_polya_anchors.parquet` with the five POD5 files. Feature metadata may come from the base extractor or its cross-boundary extension. Choose either `--inputs` or `--rebuild-inputs`.

`--reference-results R` additionally compares the 30 selected fields directly with the original `R/final_primary/host_terminal_nona_direct_scores.parquet` and `R/final_primary/mhv_reads_terminal_nona_direct.parquet`. Standard validation checks table structure, QC, exports and the expected barcode summary. `--advanced` selects the [11-stage advanced workflow](ADVANCED_ANALYSIS.md).
