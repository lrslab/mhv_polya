#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Two-stage ONT direct RNA mapping: host -> extract unmapped reads -> virus
# Input: five prepared FASTQ files
# Output per sample: host/virus sorted BAMs, BAI indices, flagstat, and summary
#
# Usage:
#   ./run_host_then_virus_5fastq_directRNA.sh HOST_FA VIRUS_FA OUTDIR fq1 fq2 fq3 fq4 fq5
#
# Example:
#   conda activate ontmap
#   ./run_host_then_virus_5fastq_directRNA.sh \
#     /s4/qqjiang/.../GRCm39/GCF_000001635.27_GRCm39_genomic.fna \
#     /s4/qqjiang/.../MHV_A59/GCF_003971785.1_ASM397178v1_genomic.fna \
#     ./aln_host_then_virus \
#     bc04.fastq bc05.fastq bc06.fastq bc07.fastq bc11.fastq
#
# Threads:
#   T_MM2=16 T_ST=8 ./run_host_then_virus_5fastq_directRNA.sh ...
###############################################################################

if [[ $# -lt 8 ]]; then
  echo "ERROR: need HOST_FA, VIRUS_FA, OUTDIR and 5 FASTQ files."
  echo "Usage: $0 HOST_FA VIRUS_FA OUTDIR fq1 fq2 fq3 fq4 fq5"
  exit 1
fi

HOST_FA="$1"
VIRUS_FA="$2"
OUTDIR="$3"
shift 3

FASTQS=("$@")
if [[ ${#FASTQS[@]} -ne 5 ]]; then
  echo "ERROR: please provide exactly 5 FASTQ files; got ${#FASTQS[@]}"
  exit 1
fi

# threads (override via env)
T_MM2="${T_MM2:-16}"
T_ST="${T_ST:-8}"

# sanity checks
command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 not found in PATH. Did you 'conda activate ontmap'?"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH. Did you 'conda activate ontmap'?"; exit 1; }

[[ -s "$HOST_FA" ]] || { echo "ERROR: HOST_FA not found or empty: $HOST_FA"; exit 1; }
[[ -s "$VIRUS_FA" ]] || { echo "ERROR: VIRUS_FA not found or empty: $VIRUS_FA"; exit 1; }

mkdir -p "$OUTDIR"/{logs,ref,host,virus,unmapped_fastq,stats}

HOST_MMI="${OUTDIR}/ref/host.dRNA.mmi"
VIRUS_MMI="${OUTDIR}/ref/virus.dRNA.mmi"

# Build indexes once
if [[ ! -s "$HOST_MMI" ]]; then
  echo "[INFO] Building host index: $HOST_MMI"
  minimap2 -d "$HOST_MMI" "$HOST_FA" > "${OUTDIR}/logs/index_host.log" 2>&1
fi

if [[ ! -s "$VIRUS_MMI" ]]; then
  echo "[INFO] Building virus index: $VIRUS_MMI"
  minimap2 -d "$VIRUS_MMI" "$VIRUS_FA" > "${OUTDIR}/logs/index_virus.log" 2>&1
fi

# Record versions
{
  echo "date: $(date -Is)"
  echo "minimap2: $(minimap2 --version 2>&1 | head -n 1)"
  echo "samtools: $(samtools --version 2>&1 | head -n 1)"
  echo "HOST_FA: $HOST_FA"
  echo "VIRUS_FA: $VIRUS_FA"
  echo "T_MM2: $T_MM2"
  echo "T_ST: $T_ST"
} > "${OUTDIR}/stats/run_info.txt"

# Run all samples
for fq in "${FASTQS[@]}"; do
  [[ -s "$fq" ]] || { echo "ERROR: FASTQ not found or empty: $fq"; exit 1; }

  base="$(basename "$fq")"
  sample="${base%%.*}"

  echo "[INFO] ===== Sample: $sample ====="

  # ---------- Step 1: map to host ----------
  host_bam="${OUTDIR}/host/${sample}.host.sorted.bam"
  host_flagstat="${OUTDIR}/stats/${sample}.host.flagstat.txt"
  host_idxstats="${OUTDIR}/stats/${sample}.host.idxstats.txt"

  minimap2 -ax splice -uf -k14 -t "$T_MM2" "$HOST_MMI" "$fq" \
    2> "${OUTDIR}/logs/${sample}.host.minimap2.log" \
    | samtools view -@ "$T_ST" -b 2> "${OUTDIR}/logs/${sample}.host.samtools_view.log" \
    | samtools sort -@ "$T_ST" -o "$host_bam" 2> "${OUTDIR}/logs/${sample}.host.samtools_sort.log"

  samtools index -@ "$T_ST" "$host_bam" 2> "${OUTDIR}/logs/${sample}.host.samtools_index.log"
  samtools flagstat -@ "$T_ST" "$host_bam" > "$host_flagstat"
  samtools idxstats "$host_bam" > "$host_idxstats"

  # ---------- Step 2: extract unmapped reads from host BAM -> FASTQ ----------
  unmapped_fq="${OUTDIR}/unmapped_fastq/${sample}.unmapped_from_host.fastq"
  samtools view -@ "$T_ST" -f 4 -b "$host_bam" \
    | samtools fastq -@ "$T_ST" - \
    > "$unmapped_fq"

  # ---------- Step 3: map unmapped reads to virus ----------
  virus_bam="${OUTDIR}/virus/${sample}.virus.sorted.bam"
  virus_flagstat="${OUTDIR}/stats/${sample}.virus.flagstat.txt"
  virus_idxstats="${OUTDIR}/stats/${sample}.virus.idxstats.txt"

  minimap2 -ax splice -uf -k14 -t "$T_MM2" "$VIRUS_MMI" "$unmapped_fq" \
    2> "${OUTDIR}/logs/${sample}.virus.minimap2.log" \
    | samtools view -@ "$T_ST" -b 2> "${OUTDIR}/logs/${sample}.virus.samtools_view.log" \
    | samtools sort -@ "$T_ST" -o "$virus_bam" 2> "${OUTDIR}/logs/${sample}.virus.samtools_sort.log"

  samtools index -@ "$T_ST" "$virus_bam" 2> "${OUTDIR}/logs/${sample}.virus.samtools_index.log"
  samtools flagstat -@ "$T_ST" "$virus_bam" > "$virus_flagstat"
  samtools idxstats "$virus_bam" > "$virus_idxstats"

  # ---------- Step 4: summary ----------
  host_total=$(awk 'NR==1{print $1}' "$host_flagstat")
  host_mapped=$(awk '/ mapped \(/ {print $1; exit}' "$host_flagstat")
  virus_total=$(awk 'NR==1{print $1}' "$virus_flagstat")
  virus_mapped=$(awk '/ mapped \(/ {print $1; exit}' "$virus_flagstat")

  # count reads in unmapped fastq (NR/4)
  unmapped_n=$(awk 'END{printf "%.0f\n", NR/4}' "$unmapped_fq" 2>/dev/null || echo "NA")

  {
    echo -e "sample\t${sample}"
    echo -e "input_fastq\t${fq}"
    echo -e "host_total_reads_in_bam\t${host_total}"
    echo -e "host_mapped_reads\t${host_mapped}"
    echo -e "unmapped_to_host_fastq_reads\t${unmapped_n}"
    echo -e "virus_total_reads_in_bam\t${virus_total}"
    echo -e "virus_mapped_reads\t${virus_mapped}"
    echo -e "host_bam\t${host_bam}"
    echo -e "virus_bam\t${virus_bam}"
    echo -e "unmapped_fastq\t${unmapped_fq}"
  } > "${OUTDIR}/stats/${sample}.summary.tsv"

done

echo "[INFO] All done."
echo "[INFO] Host BAMs : ${OUTDIR}/host/"
echo "[INFO] Virus BAMs: ${OUTDIR}/virus/"
echo "[INFO] Unmapped FASTQ: ${OUTDIR}/unmapped_fastq/"
echo "[INFO] Stats: ${OUTDIR}/stats/"
