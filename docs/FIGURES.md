# Advanced analysis figures

The plotting scripts generate Figures 1–4 from a completed `--advanced` run.
The primary barcode figure is described in [Figure 6A](FIGURE6A_METHODS.md). Each figure includes a table of its plotted values.

Each figure is supplied as:

- `png`: 400-dpi raster image for rapid review;
- `pdf`: vector output for manuscript submission;
- `svg`: editable vector output;
- `source_data/*.tsv`: the values plotted in each panel.

Recreate all figures from the package root with:

```bash
conda activate mhv-terminal-nona
python scripts/analysis/plot_publication_figures.py --package-root /path/to/advanced_analysis_run
```

The endpoint throughout is **terminal-junction non-A-like signal**, defined by
the two-sided 5 pA cutoff in [Methods](METHODS.md).

## Figure 1: per-sample current distributions

Files: `figure1_sample_current_distributions.{png,pdf,svg}`.

Panels A-E show normalized 0.25-pA histograms of the signed junction-current
deviation from the pooled long-host reference. The bc04 mock panel contains
host reads only. Panels for bc05, bc06, bc07 and bc11 overlay host-mapped and
MHV-mapped reads. The solid line marks zero after host-reference
centering, while dashed lines mark the -5 and +5 pA classification boundaries.
Panel F explains the physical and analytical scale of the cutoff.

The four infected MHV distributions have medians of -0.37, -0.66, +1.74 and
-0.50 pA in bc05, bc06, bc07 and bc11, respectively. Their interquartile widths
range from 4.64 to 5.22 pA, so the 5-pA classification distance is approximately
one full interquartile width of the observed per-library distributions. Five
picoamperes equals 5 x 10^-12 A. The rule applies to the centred contrast
between two 31-sample window means.

Source data:

- `source_data/figure1_current_distribution.tsv`
- `source_data/figure1_current_distribution_summary.tsv`

## Figure 2: high-confidence sgRNA results

Files: `figure2_sgrna_wt_mut.{png,pdf,svg}`.

Panel A shows callable high-confidence ORF4, ORF5/E, M and N sgRNA read depth
within every infected barcode. Panel B shows pooled WT and Mut read fractions,
descriptive Wilson intervals and the individual barcode points. Panel C shows
the corresponding descriptive Mut-minus-WT percentage-point difference. The
plotted group summaries are ORF4, 39/330 (11.82%) versus 73/428 (17.06%);
ORF5/E, 13/150 (8.67%) versus 29/214 (13.55%); M, 75/575 (13.04%) versus
114/813 (14.02%); and N, 344/2,788 (12.34%) versus 599/3,987 (15.02%).

All four descriptive contrasts are positive, with the largest differences for
ORF4 (+5.24 percentage points) and ORF5/E (+4.88 points). sgRNA2 and S are
omitted from the main comparison because each genotype has at most 10 callable
reads for those transcripts. The panels report descriptive pooled estimates,
read-level intervals and values from the two barcodes per genotype.

Source data:

- `source_data/figure2_sgrna_by_barcode.tsv`
- `source_data/figure2_sgrna_by_group.tsv`

## Figure 3: score behavior and positive-read lengths

Files: `figure3_signal_and_positive_read_lengths.{png,pdf,svg}`.

Panel A is the empirical survival curve of the absolute deviation from the
pooled host A|DNA-junction reference for callable MHV-mapped reads from the four
infected libraries. The vertical line is the primary 5-pA cutoff; its
intersections reproduce the 13.65% WT and 16.12% Mut read-pooled results.
Panels B and C give the median and interquartile range of poly(A)-tail length
and read length, respectively, for positive high-confidence ORF4, ORF5/E, M
and N sgRNA reads.

The read-length medians follow the expected transcript ordering. Among the
plotted positive reads, 96.0% of WT and 54.2% of Mut poly(A) lengths come from
Dorado `pt`; the remaining values use the WDX rate-scaled estimate. The source
mix is provided alongside the length summaries.

Source data:

- `source_data/figure3_signal_survival.tsv`
- `source_data/figure3_positive_length_summary.tsv`
- `source_data/figure3_polya_length_source_mix.tsv`

## Figure 4: polarity and cutoff sensitivity

Files: `figure4_direction_and_cutoff_sensitivity.{png,pdf,svg}`.

Panel A decomposes the primary two-sided call into negative- and
positive-direction deviations for mock host reads and for host and MHV reads
in each infected barcode. The same barcode-associated polarity occurs in host
and MHV reads. For example, bc06 is mainly negative-direction (host 14.45%, MHV
12.33%) whereas bc07 is mainly positive-direction (host 19.01%, MHV 15.35%).
The panel describes the barcode-level direction of the two-sided endpoint.

Panel B shows descriptive 3-, 5- and 7-pA sensitivity results; the thin lines
are individual barcodes and the thick lines are read-pooled groups. Panel C
shows that the pooled Mut-minus-WT difference changes from +5.00 percentage
points at 3 pA to +2.48 at the primary 5-pA cutoff and +0.16 at 7 pA. The 3-
and 7-pA results describe sensitivity around the primary 5-pA cutoff.

Source data:

- `source_data/figure4_signed_polarity.tsv`
- `source_data/figure4_cutoff_sensitivity_by_barcode.tsv`
- `source_data/figure4_cutoff_sensitivity_pooled.tsv`
- `source_data/figure4_cutoff_difference.tsv`

## Interpretation

The figures describe pooled read estimates alongside barcode-level variation.
The endpoint and statistical interpretation are defined in [Methods](METHODS.md).
