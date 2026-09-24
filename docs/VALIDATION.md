# Result validation

The figure workflow ends with `validate_figure.py`. It checks:

- Host and MHV table dimensions, unique read IDs, missing calls and the inclusive 5 pA cutoff.
- Agreement between exported read tables and internal scores.
- The four-barcode summary against the expected study counts and fractions.
- The pooled host reference: 26,997 reads and −10.190704345703125 pA.
- PNG, PDF and SVG outputs.

These form ten checks. Expected counts are in [the barcode summary](../expected/figure6A/figure6A_overall_virus_barcode_summary.tsv); table definitions are in [figure_io_schema.json](../schemas/figure_io_schema.json).

## Repeat the checks

From the repository directory:

```bash
python scripts/validate_figure.py --run /path/to/barcode_figure_run
python scripts/check_bundle.py
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s scripts -p 'test_*.py' -v
```

To compare all 30 retained fields with archived per-read results, supply their results directory:

```bash
python scripts/validate_figure.py \
  --run /path/to/barcode_figure_run \
  --reference-results /path/to/publication_terminal_nona_results_v1/results
```

The test suite covers signal windows, callability, the inclusive two-sided cutoff, host-reference selection, computed figure labels, result comparisons and workflow execution.

## Recorded reconstruction

On September 23, 2026, the six-stage primary workflow completed from the 36 declared inputs in the study server's Linux `py311` environment. All retained fields matched for 82,709 host reads and 86,406 MHV reads. The PNG matched the archived server rendering pixel for pixel and was visually compared with the manuscript figure.

The optional advanced reconstruction completed all 11 stages and 40 result checks. Its full 74-column host table, 82-column MHV table, four compressed exports and 34 summaries matched the study outputs.

Local checks on September 24 exercised the validators against archived outputs and passed the 15-test suite. They verify the current code and table comparisons; the complete source-input run is the server reconstruction above.

## Publication version

Complete the data accession in the README, author/version/repository information in `CITATION.cff`, and the selected `LICENSE`. Validate the full workflow in a newly installed recipe environment before archiving the publication release. [Environment details](REPRODUCIBILITY.md) and [GitHub/Zenodo instructions](GITHUB_ZENODO.md) describe these steps.
