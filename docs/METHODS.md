# Methods

## Direct RNA nanopore sequencing and sample groups

Direct RNA libraries were sequenced on an Oxford Nanopore Technologies MinION Mk1B using an RNA004 flow cell (FLO-MIN004RA) and the SQK-RNA004 chemistry. A DNA splint/adaptor carrying a sample-specific barcode was ligated at the RNA 3′ end. The five analysed barcode groups were bc04 (mock), bc05 and bc06 (MHV wild type), and bc07 and bc11 (MHV nsp15 mutant). Raw ionic-current data were stored in POD5 format and sampled at 4,000 Hz. Reads were basecalled with Dorado v1.4.0 using `rna004_130bps_sup@v5.3.0` and `--estimate-poly-a`.

## Read mapping and definition of the analysis universe

Basecalled reads were aligned sequentially, first to the *Mus musculus* GRCm39 reference and then, for reads unmapped to the host, to the MHV-A59 genome (NC_048217.1; 31,335 nt). Minimap2 v2.30-r1287 was run with the direct-RNA splice-aware settings `-ax splice -uf -k14`; alignments were coordinate-sorted and indexed with SAMtools v1.23. The metadata join used primary alignments to define host/MHV membership; its barcode and mapping rules are specified in [INPUT_OUTPUT.md](INPUT_OUTPUT.md).

Adaptor and poly(A)-boundary coordinates were obtained with WarpDemuX using the WDX6 RNA004 barcode configuration (`WDX6_rna004_v1_0`). The analysis universe was the intersection of primary MHV mappings and reads with a successful WarpDemuX boundary record: 86,406 unique MHV-mapped reads (bc04, 1,506; bc05, 11,405; bc06, 12,182; bc07, 17,284; and bc11, 44,029). The MHV BAM files contained 89,177 primary mappings before this intersection. Host-mapped reads with successful boundary records (n = 82,709) were retained to define the empirical host junction reference.

## Fixed-window terminal-junction signal statistic

The terminal endpoint was calculated from calibrated POD5 current in picoamperes (`signal_pa`). It compares the terminal RNA–DNA junction with a downstream stable-poly(A) region of the same read, using an empirical host reference to centre the contrast.

For read *i*, let \(W_i\) denote the WarpDemuX `polya_start` coordinate in raw-signal samples. A common RNA dwell approximation was fixed from the acquisition rate and nominal RNA004 translocation speed:

\[
d = 4000/130 = 30.7692307692\ \mathrm{samples\ nt^{-1}}.
\]

The integer window width was \(w=\mathrm{round}(d)=31\) samples, and the fixed half-nucleotide phase offset was \(o=\mathrm{round}(0.5d)=15\) samples. Thus the right coordinate was \(R_i=W_i+15\). The terminal-junction and stable-poly(A) means were

\[
T_i=\operatorname{mean}\{I_i[W_i-16:W_i+15)\}
\]

and

\[
A_i=\operatorname{mean}\{I_i[W_i+46:W_i+77)\},
\]

respectively, where intervals use zero-based, half-open sample coordinates. The intervening 31-sample dwell was skipped to reduce contamination of the stable-A window by the boundary transition. The signed within-read contrast was

\[
x_i=T_i-A_i.
\]

Window means and within-read contrasts were calculated using float32 calibrated signal. The same coordinate rule and dwell were used for every barcode. WDX coordinates define the signal windows; Dorado tags provide diagnostic metadata and poly(A)-length estimates.

## Empirical host A|DNA-junction reference and fixed classification threshold

An empirical A|DNA-junction reference was constructed from callable host-mapped reads with a signal-clipped WDX poly(A) span (`actual_polya_len_samples`) equivalent to at least 100 nt under the fixed dwell and a Dorado mean read quality score of at least 12. Pooling all five barcodes yielded 26,997 reference reads. The reference centre was the pooled median,

\[
x_{A}=\operatorname{median}(x_i)=-10.1907043457\ \mathrm{pA}.
\]

For every callable read, the terminal-junction deviation score was

\[
D_i=|x_i-x_A|.
\]

Reads with \(D_i\geq5.0\) pA were classified as `terminal_junction_nonA_like`; reads with \(D_i<5.0\) pA were classified as `A_junction_like`. The absolute deviation combines both current directions. The primary cutoff was 5 pA.

## Signal callability and per-read output

The endpoint used primary CNN-derived WDX boundaries. A read was callable when (i) its calibrated POD5 signal was present; (ii) mapping and WDX barcodes agreed; (iii) its boundary source was CNN; (iv) both fixed windows contained finite signal; and (v) `polya_end-polya_start` was at least 77 samples. Reads failing a criterion were retained as `uncallable`, with a missing binary call and a recorded QC reason. MHV-read callability across infected barcodes ranged from 80.44% to 83.06%.

## Descriptive summaries and reproducibility

The figure reports each infected barcode separately. Non-A-like fractions use callable MHV reads as the denominator; uncallable reads remain in the source table with a missing call. The four bars contain 1,201/9,370 (bc05), 1,415/9,799 (bc06), 2,441/14,061 (bc07) and 5,722/36,569 (bc11) positive/callable reads. Barcode labels and percentages are generated from the same summary that supplies the bars.

The reconstruction used Python 3.11.15, NumPy 2.4.4, pandas 3.0.2, pyarrow 22.0.0, POD5 0.3.44, pysam 0.24.0, SciPy 1.17.1 and Matplotlib 3.11.1. [INPUT_OUTPUT.md](INPUT_OUTPUT.md), [REPRODUCIBILITY.md](REPRODUCIBILITY.md) and [VALIDATION.md](VALIDATION.md) provide the file specifications, processing records and validation results.

## Interpretation

Non-A-like denotes an empirical junction-current deviation. Nucleotide-specific interpretation requires matched terminal-sequence standards. WT/Mut comparisons are descriptive because each genotype has two barcodes and genotype is confounded with barcode identity. Read-level intervals describe sampling among reads; biological inference depends on independent biological replication.

## Advanced analyses

Transcript assignment, tail/read-length summaries and threshold sensitivity are documented separately in [Advanced analysis](ADVANCED_ANALYSIS.md).
