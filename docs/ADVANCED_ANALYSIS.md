# Advanced analysis

This optional workflow retains the broader study analyses: canonical sgRNA and conservative gRNA assignments, host/MHV summaries, tail/read-length comparisons, signed-deviation and cutoff sensitivity, Figures 1–4, and pooled WT/Mut Figure 6A. The primary barcode figure uses the same signal scorer and does not depend on transcript assignments.

## Run

```bash
bash scripts/run_pipeline.sh \
  --project-root /path/to/MHV_polyA \
  --output /path/to/advanced_analysis_run --advanced --execute
```

The 41 source files are the primary workflow's 36 inputs plus five viral BAM indices (`bcXX.virus.sorted.bam.bai`). Use a new output directory. [Advanced input/output specification](ADVANCED_INPUT_OUTPUT.md) defines the complete file and column inventory, and [advanced_inputs.tsv](../config/advanced_inputs.tsv) lists input paths and study file sizes.

## Workflow and outputs

The 11 stages join mapping/WDX metadata, extract base and cross-boundary features, extract Dorado tags, assign transcripts, calculate the shared signal endpoint, summarize host/virus/sgRNA results, check study counts, draw Figures 1–4, draw both Figure 6A variants, and compare outputs with the frozen results.

Results are written to `results/final_primary/` and `results/final_host_virus_sgrna/`. Figures use `figures/` and `figure6A_overall_mut_vs_wt/`. Scientific exports comprise 34 TSV summaries/source tables, four compressed read tables and 18 figure files. [Figure descriptions](FIGURES.md) identify each panel.

The shared calculation is in `scripts/analysis/call_terminal_nona_direct_cutoff.py::score_all_reads`. Primary output uses `score_terminal_current.py`; the advanced caller then adds transcript annotations. Transcript assignment is implemented by `scripts/upstream/build_sgrna_master.py`, and strict sgRNA summaries by `scripts/analysis/summarize_host_virus_sgrna_direct.py`.

`sgrna_primary_high_confidence` covers canonical sgRNA and conservative gRNA. The `final_host_virus_sgrna/virus_high_conf_sgrna_*` outputs contain strict canonical sgRNA results. Signal callability and transcript identifiability are separate properties.

## MHV genomic and subgenomic RNA assignment

Diagnostic leader–body junctions were identified from CIGAR `N` operations using zero-based, half-open reference coordinates. For each `N`, the donor was the reference cursor immediately before the skipped region and the acceptor was donor plus the `N` length. Leader-side candidates required `N` ≥1,000 nt and donor ≤100. If several candidates were present, the junction nearest an expected acceptor was selected, with the longer skipped region used to break ties.

A junction was classified as canonical when MAPQ was ≥20, its acceptor was within ±30 nt of an empirical MHV alignment breakpoint—21,746 for sgRNA2/ORF2a-HE, 23,916 for sgRNA3/S, 27,931 for sgRNA4/ORF4, 28,305 for sgRNA5/ORF5-E, 28,943 for sgRNA6/M, and 29,648 for sgRNA7/N—and at least 15 aligned reference bases occurred in each 30-nt flank, `[donor-30, donor)` and `[acceptor, acceptor+30)`. The acceptor coordinates are the dominant observed alignment breakpoints.

Genomic RNA was assigned only to reads with no CIGAR `N`, MAPQ ≥20, and at least 100 aligned bases in the ORF1-unique interval `[210, 21743)`. Reads failing both the canonical sgRNA and conservative gRNA criteria were retained as `nested_RNA_unidentifiable`, because common 3′ coverage from 5′-truncated direct-RNA reads cannot identify the nested transcript of origin. Nearest-body labels were retained only as low-confidence sensitivity annotations and were excluded from primary sgRNA summaries. ORF1-positive reads containing a noncanonical skipped region were flagged as possible defective viral genome or chimeric reads rather than assigned as gRNA.

The high-confidence transcript set comprised 12,617 canonical sgRNA reads and 52 conservative gRNA reads. The canonical sgRNA counts were 3 sgRNA2/ORF2a-HE, 30 S, 1,010 ORF4, 498 ORF5/E, 1,867 M, and 9,209 N reads. The remaining 73,737 MHV reads were transcript-unidentifiable. gRNA and unassigned reads have separate summaries.

## Poly(A)-tail and read lengths

Poly(A)-tail length was taken preferentially from a positive Dorado `pt` BAM tag (65,481 reads). If `pt` was absent or non-positive, a WarpDemuX rate-scaled estimate was used (20,925 reads):

\[
\mathrm{poly(A)}_{nt}=\mathrm{poly(A)}_{samples}\times
\frac{\mathrm{basecalled\ read\ length}_{nt}}{\mathrm{total\ signal\ length}_{samples}}.
\]

Basecalled read length was obtained from the merged FASTQ length table and agreed exactly with BAM query length for all 86,406 reads. Poly(A)-tail and read lengths were retained for every strict-sgRNA `terminal_junction_nonA_like` read and summarized by sample, barcode, sgRNA, and signal class using medians and interquartile ranges.

## Additional signal summaries

The advanced workflow also reports 3 and 7 pA sensitivity thresholds, host-reference and MHV current distributions, signed-deviation polarity, and pooled genotype estimates. Poly(A)-length comparisons depend on the mixture of Dorado and WDX estimates, whose proportions differ between groups. Read-level intervals describe sampling among reads; genotype inference depends on biological replication.

## Validation

The 40 checks cover the structure and QC of the full 74-column host and 82-column MHV tables, four compressed read exports and 34 expected TSV comparisons. `--reference-results` adds direct numerical comparisons with archived read-level tables. Repeat the comparison with:

```bash
python scripts/validate_target.py --run /path/to/advanced_analysis_run
```

[Validation](VALIDATION.md) summarizes the recorded reconstruction and result checks.
