# Environment and upstream processing

The [README](../README.md) lists every command needed to reproduce the barcode figure. [run_figure.sh](../scripts/run_figure.sh) executes those five analysis steps and the final validation. The 36 source inputs are listed in [inputs.tsv](../config/inputs.tsv).

## Environment

| Component | Recorded version |
|---|---|
| Python | 3.11.15 |
| NumPy / pandas | 2.4.4 / 3.0.2 |
| pyarrow / POD5 / pysam | 22.0.0 / 0.3.44 / 0.24.0 |
| SciPy / Matplotlib | 1.17.1 / 3.11.1 |
| SAMtools used to read archived BAMs | 1.17 |
| SAMtools in `environment.yml` | 1.23 |

The installation recipe pins the Python packages and uses SAMtools 1.23 to resolve the environment's zlib dependency. Conda resolution and pip Linux wheel checks passed. Full dataset reconstruction used the original Linux `py311` environment; a full run in a newly installed recipe environment remains to be completed.

The inputs occupy 6,222,876,270 bytes. The primary reconstruction took 2 minutes 36 seconds and peaked at 2.54 GiB resident memory on the study server. Allow 4 GB RAM and 2 GB output space.

## Run on a server

```bash
mhv_output=/path/to/barcode_figure_run
mkdir -p "$(dirname "$mhv_output")"
nohup bash scripts/run_figure.sh /path/to/MHV_polyA "$mhv_output" \
  > "${mhv_output}.log" 2>&1 < /dev/null &
```

The workflow stores generated features under `inputs/`, intermediate tables under `work/`, and step logs under `logs/`, all inside the output directory. The scoring step reads [analysis_config.json](../config/analysis_config.json); the validator reads the primary table schema and expected barcode summary.

For cached metadata, command previews or an execution manifest, the existing `scripts/run_pipeline.sh --project-root P --output O` runner remains available. Its `--plan`, `--inputs`, `--reference-results` and `--execute` options are described in [Input/output definitions](INPUT_OUTPUT.md). The [advanced workflow](ADVANCED_ANALYSIS.md) adds transcript analysis and other study figures.

## Preparation of the study inputs

The reproducible starting point is demultiplexed POD5 plus archived WDX boundaries/predictions, Dorado BAMs and alignments. These preserve the study's read assignments and original signal coordinates.

| Processing step | Study settings |
|---|---|
| Dorado | 1.4.0+ba44a013; model `rna004_130bps_sup@v5.3.0`; ordinary basecalls and `--estimate-poly-a` BAMs |
| Alignment | Minimap2 2.30-r1287 and SAMtools 1.23 |
| Host reference | Mouse GRCm39, `GCF_000001635.27` |
| Viral reference | MHV-A59, `GCF_003971785.1`; contig `NC_048217.1` |
| WDX model | `WDX6_rna004_v1_0` |
| WDX boundary run | March 18, 2026; boundary saving enabled; `core.max_obs_trace=50000`; 8 cores; output batches of 40,000; minibatches of 1,000 |

Initial WDX demultiplexing on March 9 and boundary extraction on March 18 were separate runs. This analysis uses the March 18 boundary and prediction exports. The exact WDX source commit and complete original POD5 subset/merge command are unavailable in the retained records, so the specified preprocessing outputs are part of the deposited input set.

[run_host_then_virus_5fastq_directRNA.sh](../scripts/upstream/run_host_then_virus_5fastq_directRNA.sh) preserves the upstream alignment procedure. It builds MMI indices with default settings and aligns with `-ax splice -uf -k14`; effective indexing parameters follow the saved MMI indices.
