# Advanced analysis: input and output specification

This specification describes the optional advanced analyses using the within-read 5 pA workflow for `publication_terminal_nona_results_v1`. Run `bash scripts/run_pipeline.sh --project-root P --output O --advanced --execute`, using absolute paths for the input data directory `P` and a new output directory `O`.

The workflow reads 41 source files and produces 34 TSV summaries/source tables, four gzip-compressed read tables and 18 figure files. Parquet and NPY intermediates are generated under `O`.

The [input inventory](../config/advanced_inputs.tsv) lists the 41 source paths and study file sizes. The output tables are defined below. The [schema inventory](../schemas/observed_io_schema.json) records ordered headers, table dimensions, storage types, missing-value counts and array shapes observed in the study.

## Command-line inputs

| Argument | Definition |
|---|---|
| `--project-root P` | Required. Root containing `MHV_Multiple_sample/` with the study inputs listed below. |
| `--output O` | Required. A new directory for results, intermediates, logs and validation reports. |
| `--advanced` | Select the full 11-stage workflow with transcript analysis. |
| `--rebuild-inputs` | Explicitly select the default reconstruction from 41 source files. |
| `--plan` | Preview the stage commands and input paths. WDX paths use wildcards when the input directory is unavailable. |
| Neither `--plan` nor `--execute` | Check file presence, dependencies and output-path availability. |
| `--execute` | Execute after preflight; stop at the first failed stage with a nonzero exit status. |
| `--inputs C` | Optional internal cache mode: use 11 existing compact Parquet tables plus the five POD5 files. Mutually exclusive with `--rebuild-inputs`; see below. |
| `--reference-results R` | Optional direct comparison with archived reference files, in addition to the bundled target checks. |

`config/analysis_config.json` records the parameters implemented by the runner and caller. Use `--rebuild-inputs` to start from source data or `--inputs` to select a compact cache. Omitting both rebuilds intermediates from the source files.

## Direct analysis inputs: 41 files

All paths below are relative to `P`. Expand `bcXX` to `bc04`, `bc05`, `bc06`, `bc07`, and `bc11`. Expand `{0..6}` to the seven WDX shard numbers, 0 through 6.

Define these path prefixes:

```text
S = MHV_Multiple_sample
D = S/20260306_5bar/20260306_1823_MN25294_FAZ51928_089bb3a7/pod5_by_barcode
B = S/Bam/aln_host_then_virus
W = S/Warpdemux_polyA/warpdemux_WDX6_rna004_v1_0_20260318_1039_ade68d6a
```

| Files | Count | Format and purpose |
|---|---:|---|
| `D/bcXX.pod5` | 5 | Demultiplexed raw signal and calibration; supplies per-read current in pA. |
| `B/host/bcXX.host.sorted.bam` | 5 | Host-first primary alignments; host membership, barcode, reference and MAPQ. |
| `B/virus/bcXX.virus.sorted.bam` | 5 | Viral primary alignments; membership and CIGAR-based transcript assignment. |
| `B/virus/bcXX.virus.sorted.bam.bai` | 5 | Matching indices required for indexed access to the viral BAMs. |
| `S/polyA_bam/polya_5pod5_results/bcXX.polya.bam` | 5 | Archived Dorado poly(A) estimates, quality scores and diagnostic signal tags. |
| `W/boundaries/detected_boundaries_{0..6}.csv` | 7 | WDX boundaries in original raw-signal sample coordinates. |
| `W/predictions/barcode_predictions_{0..6}.csv.gz` | 7 | WDX barcode predictions from the same archived March 18 run. |
| `S/Bam/fastq.read_lengths.merged.tsv` | 1 | Basecalled read lengths in nucleotides. |
| `W/boundaries/polya_nt_estimated.tsv` | 1 | Archived WDX length estimates used for metadata QA. |
| **Total** | **41** | Source inputs for the study reconstruction. |

The input inventory records study paths and file sizes. The runner discovers WDX shards by filename pattern; the study uses seven boundary and seven prediction shards. Preflight checks file presence, and final validation compares the numerical results. Demultiplexing, basecalling and alignment precede this workflow; their records are summarized in [Reproducibility](REPRODUCIBILITY.md).

Optional `W/failed_reads/failed_reads_*.csv.gz` exports contain `read_id` and optionally `fail_reason`. They contribute failed-read QA counts under `O/work/`. The study reconstruction produced the same analysis metadata from the 41 declared source inputs.

### Shared identifiers, groups and coordinates

Every join uses the original read identifier: POD5 `read_id`, BAM query name and text-table `read_id` refer to the same UUID. Mapping barcodes are `bc04`, `bc05`, `bc06`, `bc07`, `bc11`; group labels are `mock`, `MHV-wt`, `MHV-wt`, `MHV-mut`, `MHV-mut`, respectively. WDX numeric barcode 4 is converted to `bc04`, and similarly for the other samples.

Signal coordinates are zero-based sample indices into the original POD5 signal, with half-open windows. The recorded experiment is RNA004 at 4,000 samples/s. WDX coordinates must correspond to that signal before basecaller trimming. BAM/reference coordinates used for transcript assignment are zero-based, half-open after conversion by pysam. The viral reference contig is `NC_048217.1`.

### POD5 and alignment requirements

POD5 supplies calibration and original read UUIDs. The caller reads `ReadRecord.signal_pa` as float32 calibrated pA for the selected reads. A missing selected signal is recorded as uncallable.

The mapping join reads BAM records with `samtools view -F 2304`, excluding secondary and supplementary alignments. It reads query name, flag, reference name and MAPQ. Unmapped primary records contribute barcode membership but not mapped host/virus status. If several mapped candidates exist, the join first prefers the WDX-matching barcode, then host over virus, then larger MAPQ. This preserves the historical host-first mapping rule.

The indexed viral BAMs provide CIGAR, query length, reference position and MAPQ for transcript assignment. Each selected viral read has one primary alignment, and each BAI indexes its corresponding BAM. [Advanced analysis](ADVANCED_ANALYSIS.md) defines the six canonical sgRNA acceptors and conservative gRNA rule.

Dorado poly(A) BAMs are scanned sequentially. Read IDs must be unique within and across the five files. Tags have the following roles:

| Tag | Meaning in this workflow | Missing-value behavior |
|---|---|---|
| `qs` | Dorado mean read quality score, used for host reference selection | Missing becomes null and excludes a host read from the reference. |
| `pt` | Poly(A)-length estimate, nt | A positive value is preferred; missing/nonpositive values use the WDX estimate. |
| `pa` | Array of search anchor, primary start/end, and optional secondary start/end in signal samples | Diagnostic metadata; missing/short arrays are padded with nulls. |
| `ts`, `ns` | Trimmed signal start and basecalled signal end | Optional, retained as nullable metadata. |
| `sm`, `sd` | Dorado scaling midpoint and dispersion | Optional, retained as nullable metadata. |

### WDX boundary CSV schema

The original CSV files have 55 columns. The scripts consume the following 19 names across boundary and anchor extraction:

```text
read_id, signal_len, adapter_start, adapter_end,
polya_start, polya_end, polya_len,
polya_mean, polya_std, polya_med, polya_mad,
llr_adapter_end, llr_polya_end,
cnn_adapter_end, cnn_polya_end,
start_peak_adapter_end, start_peak_polya_end,
adapter_dt_med, adapter_dt_mad
```

`read_id` is a string. Signal lengths and boundary positions are numeric sample counts; selected boundaries use finite integral positions. Dwell median/MAD are samples, and alternative detector positions can be missing. The caller assigns `boundary_method = cnn` when both `cnn_adapter_end` and `cnn_polya_end` are present, and `fallback` otherwise. The endpoint anchor is the exported `polya_start`. WDX `polya_mean/std/med/mad` values are retained as upstream metadata; the endpoint's calibrated window means are calculated from POD5.

Read IDs are unique across boundary shards. The join records duplicate counts and the anchor extractor checks uniqueness. The anchor table's `warpdemux_boundary_method` is a diagnostic field based on CNN/LLR adapter-end availability; the caller uses `boundary_method` for eligibility.

### Prediction and length-table schemas

| Input | Required or consumed columns | Interpretation |
|---|---|---|
| `barcode_predictions_*.csv.gz` | `#read_id` (or `read_id`), `predicted_barcode`; `confidence_score` when present | Comma-delimited gzip. Identifier string, numeric barcode label, numeric confidence. Negative/unclassified predictions do not identify a barcode. Additional class probabilities are unused. Duplicate predictions retain the larger confidence in the historical join. |
| `fastq.read_lengths.merged.tsv` | `sample`, `read_id`, `read_nt_len` | Tab-delimited, 226,920 data rows. Barcode, read identifier, length in nt. Selected viral lengths are nonmissing and agree with BAM query length in this study. The join retains the first duplicate ID. |
| `polya_nt_estimated.tsv` | `read_id`, `polya_nt_est` | Tab-delimited, 268,378 data rows. Numeric estimate in nt for QA comparison only. The original eight-column header is recorded in the schema inventory. |

The mapping join calculates `warpdemux_polya_estimate_nt = polya_len * read_nt_len / signal_len`, serializes it to six decimal places and uses it as the fallback poly(A) length. Missing optional detector values and tags remain null.

## Analysis units and missing-value rules

The selected host/virus metadata contains 169,115 reads: 82,709 host and 86,406 MHV. Each row is one read, and group summaries pool read counts across barcodes.

For WDX `polya_start = W`, the terminal mean uses `signal_pa[W-16:W+15]`, and the stable-A mean uses `signal_pa[W+46:W+77]`. Both windows contain 31 samples. The within-read contrast is `x = terminal - stableA`. The pooled host reference is the median x among callable host reads with Dorado `qs >= 12` and `actual_polya_len_samples >= 100 * 30.7692307692`. The latter field is the signal-clipped span, `max(0, min(polya_end, len(signal_pa)) - polya_start)`.

The reference has 26,997 reads and `x_A = -10.190704345703125` pA. The signed score is `x - x_A`, the absolute score is `D = abs(x - x_A)`, and the positive call is `D >= 5.0` pA. Three and seven pA are separate sensitivity summaries.

Callability requires mapped host/virus membership, available raw signal, concordant mapping/WDX barcodes, the primary CNN boundary, a detected span `polya_end - polya_start >= 77`, and complete finite terminal/stable-A windows. `eligible_terminal_nona_direct` records this combined rule. Other feature/QC columns are retained as metadata. A read can have finite window means while failing another eligibility criterion.

## Scientific outputs

All paths here are relative to `O`. TSVs use UTF-8, tabs and one header row, with the columns listed in the schema inventory. Empty cells represent missing values; booleans serialize as `True`/`False`. Read-level gzip tables contain the same text format after decompression. A binary call has three states: positive, negative and uncallable (missing).

### Read tables

| Path | Rows × columns | Row population |
|---|---:|---|
| `results/final_primary/mhv_reads_terminal_nona_direct.tsv.gz` | 86,406 × 82 | All selected MHV reads, including mock and uncallable reads. |
| `results/final_primary/mhv_terminal_junction_nona_like_reads.tsv.gz` | 10,962 × 82 | Callable MHV reads with `D >= 5` pA. |
| `results/final_host_virus_sgrna/virus_high_conf_sgrna_nonA_like_reads.tsv.gz` | 1,307 × 13 | Positive strict canonical sgRNA reads; excludes gRNA and transcript-unidentifiable reads. |
| `figure6A_overall_mut_vs_wt/data/figure6A_overall_virus_read_level.tsv.gz` | 84,900 × 19 | All infected-barcode MHV reads, including uncallable reads and all transcript classes; excludes bc04. |

The primary MHV table is the full read-level text export. Host results comprise TSV summaries and an internal 82,709 × 74 Parquet. Figure 6A's read table uses eight decimal places and its summaries use 12 significant digits. The primary MHV export retains fuller text precision, while the internal tables retain the values used for exact-value validation.

### Core data dictionary

Complete ordered headers for all exported tables are in the schema inventory. In the full MHV table, the first 49 columns contain inherited mapping, boundary, feature and Dorado metadata. The endpoint and transcript fields have these meanings:

| Field | Type / unit | Definition |
|---|---|---|
| `read_id` | String | Unique original sequencing read identifier. |
| `mapping_group` | Category | `host` or `virus`; the full MHV export contains `virus`. |
| `mapping_barcode`, `sample_group` | Categories | Mapping barcode and experimental group as defined above. |
| `barcode_mapping_consistent` | Numeric 0/1 | Mapping and WDX barcode agreement; stored as floating-point 0.0/1.0 in the full target table. |
| `boundary_method` | Category | `cnn` or `fallback`, from the WDX detector columns. |
| `terminal_nona_anchor_mode`, `terminal_nona_coordinate_system` | Strings | `warp_raw` and `raw_signal_sample`. |
| `terminal_nona_detected_right_array_index` | Samples | W, the WDX `polya_start` in raw POD5 coordinates. |
| `terminal_nona_model_right_array_index` | Samples | W + 15, the terminal window's exclusive right edge. |
| `terminal_nona_terminal_left_array_index` | Samples | W − 16, the terminal window's inclusive left edge. |
| `terminal_nona_pa0_array_index` | Samples, nullable | Inherited Dorado `pa` search-anchor tag; diagnostic only. |
| `terminal_nona_dorado_warp_delta_samples` | Samples, nullable | WDX `polya_start` minus the Dorado primary start tag; diagnostic only. |
| `terminal_nona_rna_dwell_samples_per_nt` | Samples/nt | Fixed dwell of 30.7692307692 for every read. |
| `terminal_nona_right_phase_nt` | Nominal nt | Fixed phase 0.5; rounded to 15 samples in the window calculation. |
| `terminal_window_mean_centered_pa` | pA | Absolute calibrated mean of the terminal raw window. |
| `stableA_window_mean_centered_pa` | pA | Absolute calibrated mean of the stable-A raw window. |
| `stableA_window_mad_pa` | pA | Median absolute deviation within the stable-A window. |
| `terminal_minus_stableA_pa` | pA | Within-read contrast x; can be finite on an otherwise uncallable read. |
| `eligible_terminal_nona_direct` | Integer 0/1 | Combined signal-callability indicator. |
| `terminal_nona_direct_qc_status` | Category | `eligible`, `pod5_signal_not_found`, `barcode_mapping_mismatch`, `nonprimary_boundary_method`, `warp_polya_too_short_for_stableA_window`, or `terminal_or_stableA_window_incomplete`. One reason is retained by the original QC precedence. |
| `pooled_long_host_A_DNA_reference_pa` | pA | Pooled host RNA A\|DNA-junction reference x_A, repeated on every read. |
| `terminal_nona_fixed_cutoff_pa` | pA | Primary absolute-deviation cutoff, 5.0. |
| `terminal_nona_signed_difference_pa` | pA, nullable | x − x_A for callable reads; missing otherwise. |
| `terminal_nona_abs_difference_pa` | pA, nullable | Absolute signed difference for callable reads; missing otherwise. |
| `terminal_junction_nonA_like` | Nullable Boolean | `True` if callable and D ≥ 5; `False` if callable and D < 5; missing if uncallable. |
| `terminal_nona_direct_class` | Category | `terminal_junction_nonA_like`, `A_junction_like`, or `uncallable`. |
| `long_host_reference_member` | Integer 0/1 | Membership of the 26,997-read host reference; 0 for every MHV read. |
| `sgrna_primary_assignment` | Category | `sgRNA2_ORF2a_HE`, `sgRNA3_S`, `sgRNA4_ORF4`, `sgRNA5_ORF5_E`, `sgRNA6_M`, `sgRNA7_N`, `gRNA_ORF1ab`, or `sgRNA_unassigned`. |
| `sgrna_all_read_identifiability` | Category | Transcript-identifiability class, separate from signal callability. |
| `sgrna_primary_high_confidence` | Boolean | Includes 12,617 canonical sgRNA and 52 conservative gRNA reads. |
| `sgrna_low_confidence_gene_guess` | Nullable category | Sensitivity annotation; excluded from primary strict sgRNA summaries. |
| `sgrna_assignment_method`, `sgrna_assignment_confidence` | Categories | Evidence/rule and confidence for the transcript annotation. |
| `polya_length_nt`, `polya_length_source` | nt / category | Positive Dorado `pt` preferred (`dorado_pt`, 65,481 MHV reads), otherwise recalculated WDX length (`warpdemux_rate_scaled_fallback`, 20,925 MHV reads). |
| `read_nt_len` | nt | Basecalled read length from the merged length table. |

`final_primary/mhv_sgrna_*` summaries use the combined high-confidence flag and include conservative gRNA. `final_host_virus_sgrna/virus_high_conf_sgrna_*` contains the strict sgRNA results. Signal callability and transcript identifiability are separate: a transcript-unidentifiable read can have a valid signal call.

### Summary denominators

For each reported population, `n_total` counts all rows, `n_callable` counts eligible rows, and `n_nonA_like` (or the file's equivalent named count) counts positive callable rows. The non-A-like fraction is `n_nonA_like / n_callable`; callable fraction is `n_callable / n_total`; the all-read lower bound is `n_nonA_like / n_total`. A-like count is `n_callable - n_nonA_like`. A zero denominator yields missing values.

The non-A-like fraction is named `terminal_junction_nonA_like_fraction_callable` in primary summaries, `nonA_like_fraction_callable` in split summaries, and `terminal_junction_nonA_like_fraction` in Figure 6A. Fractions are on a 0–1 scale unless a field says percent or percentage points. Wilson intervals describe read-level 95% intervals. Group estimates pool reads across barcodes. `long_host_*` summaries use the 26,997-read reference subset. Figure 6A pools bc05/bc06 as WT and bc07/bc11 as Mut.

### Figures

The following six stems each produce PDF, PNG and SVG, for 18 files:

```text
figures/figure1_sample_current_distributions
figures/figure2_sgrna_wt_mut
figures/figure3_signal_and_positive_read_lengths
figures/figure4_direction_and_cutoff_sensitivity
figure6A_overall_mut_vs_wt/output/figure6A_overall_mut_vs_wt_virus
figure6A_overall_mut_vs_wt/output/figure6A_barcode_level_mut_vs_wt_virus
```

## Internal files and run reports

`O/inputs/` contains 16 rebuilt Parquet tables and 15 float32 NPY arrays: five base-feature tables (39 columns), five cross-boundary tables (44 columns), five Dorado tables (22 columns), and one rebuilt transcript master. Each barcode's base/cross-boundary rows match; the barcode row counts are 21,550, 21,962, 25,060, 28,427 and 72,116. There are two 360-sample base arrays and one 600-sample cross-boundary array per barcode; the latter's boundary index is 240.

`O/results/final_primary/` contains the 82,709 × 74 host and 86,406 × 82 MHV Parquets used by the plotting and validation stages.

The two result directories contain five JSON files in total: primary `analysis_parameters.json`, `validation_report.json`, `result_provenance.json`, and split `validation_report.json`, `result_provenance.json`, plus a README in each. The output root additionally contains:

| File or directory | Content |
|---|---|
| `run_manifest.json` | Full command plan, input paths/sizes, package versions, stage timestamps and exit codes, overall status. |
| `analysis_validation.json` | Original study count/QC checks; optional file comparisons only when external references were supplied. |
| `local_target_validation.json` | Forty checks of read-table schemas/QC and expected TSV values; optional direct reference comparisons. |
| `logs/` | Eleven per-stage logs in rebuild mode; six in compact-input mode. |
| `work/` | Mapping-join QA tables, a project-view symlink to source data and plotting cache. |

A successful run records `status: PASS` in all three root JSON reports.

## Result validation

`schemas/observed_io_schema.json` records the expected read-table dimensions and columns. Validation checks unique read identifiers, callability and the numerical summaries in `expected/`.

The 40 checks comprise two score-table schema/QC checks, four compressed read-table schema checks and 34 expected TSV comparisons (`rtol=0`, `atol=1e-12`). Figures are regenerated from the result tables.

`--reference-results R` adds direct comparisons with archived files. `R` contains `final_primary/` with the two full Parquets, two read-level gzip exports and 12 TSV summaries; and `final_host_virus_sgrna/` with the strict-positive gzip and eight TSV summaries.

## Optional compact cache

`--inputs C` selects five `C/features/bcXX.metadata.parquet`, five `C/anchors/bcXX.dorado_polya_anchors.parquet`, and `C/sgrna/mhv_sgrna_assignment_length_master.parquet`, together with the five barcode POD5 files. The historical transcript master has 86,406 rows and 38 columns. The rebuilt master uses `O/inputs/sgrna_master.parquet` and can retain additional upstream columns; the endpoint transfers its selected metadata fields. These two input layouts use different master-table paths.

The schema inventory uses `compact_input`, `rebuild_input`, `rebuilt_intermediate` and `output` scopes to identify each table's role, columns and observed missing values.

## Complete summary-table inventory

The tables below enumerate all 34 exported TSV summaries and figure-source tables. Dimensions exclude the header row and correspond to the frozen study target. Complete column lists are available in the schema inventory.

### Primary summaries (12 files)

Directory: `results/final_primary/`.

| Filename | Rows | Columns |
|---|---:|---:|
| `host_long_polya_pooled_reference_summary.tsv` | 6 | 5 |
| `long_host_barcode_nona_direct_reference_summary.tsv` | 5 | 12 |
| `long_host_group_nona_direct_reference_summary.tsv` | 3 | 11 |
| `mhv_barcode_boundary_method_nona_direct_summary.tsv` | 10 | 13 |
| `mhv_barcode_nona_direct_summary.tsv` | 5 | 12 |
| `mhv_barcode_nona_direct_threshold_sensitivity.tsv` | 15 | 13 |
| `mhv_barcode_vs_long_host_nona_direct_comparison.tsv` | 5 | 16 |
| `mhv_group_nona_direct_summary.tsv` | 3 | 11 |
| `mhv_group_vs_long_host_nona_direct_comparison.tsv` | 3 | 15 |
| `mhv_sgrna_barcode_nona_direct_summary.tsv` | 31 | 13 |
| `mhv_sgrna_group_nona_direct_summary.tsv` | 18 | 12 |
| `mhv_sgrna_nona_direct_length_summary.tsv` | 31 | 12 |

### Host/virus and strict transcript summaries (8 files)

Directory: `results/final_host_virus_sgrna/`.

| Filename | Rows | Columns |
|---|---:|---:|
| `barcode_host_vs_virus_nona_summary.tsv` | 10 | 15 |
| `virus_high_conf_grna_by_barcode_summary.tsv` | 4 | 13 |
| `virus_high_conf_sgrna_by_barcode_summary.tsv` | 27 | 13 |
| `virus_high_conf_sgrna_by_group_summary.tsv` | 16 | 12 |
| `virus_high_conf_sgrna_length_by_barcode.tsv` | 46 | 13 |
| `virus_transcript_accounting_by_barcode.tsv` | 5 | 16 |
| `virus_transcript_identifiability_by_barcode.tsv` | 36 | 13 |
| `virus_unassigned_by_barcode_summary.tsv` | 5 | 13 |

### Figures 1–4 source data (11 files)

Directory: `figures/source_data/`.

| Filename | Rows | Columns |
|---|---:|---:|
| `figure1_current_distribution.tsv` | 1,080 | 8 |
| `figure1_current_distribution_summary.tsv` | 9 | 16 |
| `figure2_sgrna_by_barcode.tsv` | 16 | 14 |
| `figure2_sgrna_by_group.tsv` | 8 | 14 |
| `figure3_polya_length_source_mix.tsv` | 4 | 4 |
| `figure3_positive_length_summary.tsv` | 8 | 13 |
| `figure3_signal_survival.tsv` | 1,202 | 4 |
| `figure4_cutoff_difference.tsv` | 3 | 4 |
| `figure4_cutoff_sensitivity_by_barcode.tsv` | 12 | 14 |
| `figure4_cutoff_sensitivity_pooled.tsv` | 6 | 5 |
| `figure4_signed_polarity.tsv` | 9 | 12 |

### Figure 6A summaries (3 files)

Directory: `figure6A_overall_mut_vs_wt/data/`.

| Filename | Rows | Columns |
|---|---:|---:|
| `figure6A_overall_virus_barcode_summary.tsv` | 4 | 13 |
| `figure6A_overall_virus_effect_summary.tsv` | 1 | 18 |
| `figure6A_overall_virus_group_summary.tsv` | 2 | 12 |
