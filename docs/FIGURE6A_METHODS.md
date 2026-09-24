# Barcode Figure 6A

## Figure caption

**A. Barcode-level terminal-junction Non-A-like fractions in MHV RNA.** The left panel shows A-junction-like (grey) and Non-A-like (blue) reads as percentages of callable MHV reads in each barcode. The right panel expands the Non-A-like fraction. bc05 and bc06 are WT-associated libraries; bc07 and bc11 are Mut-associated libraries. Non-A-like denotes an absolute deviation of at least 5 pA from the pooled host reference after subtracting stable-poly(A) current from terminal-junction current within each read. Counts are 1,201/9,370 (12.82%), 1,415/9,799 (14.44%), 2,441/14,061 (17.36%) and 5,722/36,569 (15.65%), respectively. Each barcode uses its own callable-read denominator. The comparison is descriptive.

## Method and source data

The terminal mean uses raw-signal samples `[W−16, W+15)` and the stable-poly(A) mean uses `[W+46, W+77)`, where W is the primary CNN WDX `polya_start`. Both means use calibrated POD5 current in pA. The contrast is centred on 26,997 callable host reads pooled across all five barcodes, each with a signal-clipped poly(A) span equivalent to at least 100 nominal nt and Dorado mean quality ≥12. The reference median is −10.190704345703125 pA. [Methods](METHODS.md) gives the callability rules and numerical definition.

The figure contains all callable MHV transcript classes within each infected barcode. Its source table retains all 84,900 selected infected MHV reads, of which 69,799 are callable and 10,779 are Non-A-like. Uncallable reads have a recorded QC reason and a missing call. The four-row barcode summary supplies both the bars and text labels.

The default pipeline writes the figure in PDF, SVG and 400 dpi PNG under `figure6A/output/`, with source data under `figure6A/data/`. The PDF embeds TrueType fonts. Vector exports are suitable for resizing to manuscript dimensions; check text size at the journal's final panel width.

Pooled genotype estimates and sgRNA analyses are described in [Advanced analysis](ADVANCED_ANALYSIS.md).
