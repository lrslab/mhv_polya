#!/usr/bin/env python3
"""CIGAR-based transcript assignment functions from
analyze_mhv_sgrna_terminal_u.py from the original study analysis."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pysam

SGRNA_CENTERS = {
    "sgRNA2_ORF2a_HE": 21746,
    "sgRNA3_S": 23916,
    "sgRNA4_ORF4": 27931,
    "sgRNA5_ORF5_E": 28305,
    "sgRNA6_M": 28943,
    "sgRNA7_N": 29648,
}
ORF1_UNIQUE_INTERVAL = (210, 21743)

def overlap_length(blocks: list[tuple[int, int]], start: int, end: int) -> int:
    return int(
        sum(max(0, min(block_end, end) - max(block_start, start))
            for block_start, block_end in blocks)
    )


def introns_from_alignment(alignment: pysam.AlignedSegment) -> list[dict[str, int]]:
    cursor = int(alignment.reference_start)
    introns: list[dict[str, int]] = []
    for operation, length in alignment.cigartuples or []:
        if operation == 3:  # N: skipped reference sequence
            introns.append(
                {"donor0": cursor, "acceptor0": cursor + int(length), "length": int(length)}
            )
        if operation in {0, 2, 3, 7, 8}:  # reference-consuming CIGAR ops
            cursor += int(length)
    return introns


def nearest_sgrna(position: int) -> tuple[str, int]:
    label, center = min(
        SGRNA_CENTERS.items(), key=lambda item: abs(int(position) - item[1])
    )
    return label, int(position) - center


def first_body_position(blocks: list[tuple[int, int]]) -> float:
    positions = [max(start, 201) for start, end in blocks if end > 201]
    return float(min(positions)) if positions else float("nan")


def classify_alignment(
    alignment: pysam.AlignedSegment,
    args: argparse.Namespace,
) -> dict[str, object]:
    blocks = [(int(start), int(end)) for start, end in alignment.get_blocks()]
    introns = introns_from_alignment(alignment)
    body_pos = first_body_position(blocks)
    orf1_overlap = overlap_length(blocks, *ORF1_UNIQUE_INTERVAL)

    leader_introns = [
        intron
        for intron in introns
        if intron["length"] >= args.min_intron
        and intron["donor0"] <= args.leader_donor_max
    ]
    chosen: dict[str, int] | None = None
    chosen_label: str | None = None
    chosen_delta: int | None = None
    left_flank = 0
    right_flank = 0
    if leader_introns:
        candidates: list[tuple[int, int, dict[str, int], str, int, int]] = []
        for intron in leader_introns:
            label, delta = nearest_sgrna(intron["acceptor0"])
            left = overlap_length(
                blocks,
                intron["donor0"] - args.junction_flank_window,
                intron["donor0"],
            )
            right = overlap_length(
                blocks,
                intron["acceptor0"],
                intron["acceptor0"] + args.junction_flank_window,
            )
            candidates.append((abs(delta), -intron["length"], intron, label, left, right))
        _, _, chosen, chosen_label, left_flank, right_flank = min(
            candidates, key=lambda value: (value[0], value[1])
        )
        chosen_delta = chosen["acceptor0"] - SGRNA_CENTERS[chosen_label]

    canonical_junction = bool(
        chosen is not None
        and abs(int(chosen_delta)) <= args.acceptor_tolerance
        and left_flank >= args.min_junction_flank_aligned
        and right_flank >= args.min_junction_flank_aligned
        and alignment.mapping_quality >= args.min_mapq
    )
    grna_supported = bool(
        len(introns) == 0
        and alignment.mapping_quality >= args.min_mapq
        and orf1_overlap >= args.min_orf1_overlap
    )

    if canonical_junction:
        primary = str(chosen_label)
        best_guess = primary
        method = "canonical_leader_body_junction"
        confidence = "high"
    elif grna_supported:
        primary = "gRNA_ORF1ab"
        best_guess = primary
        method = "unique_ORF1_alignment_without_N"
        confidence = "high"
    else:
        primary = "sgRNA_unassigned"
        if orf1_overlap >= args.min_orf1_overlap:
            best_guess = "ORF1_positive_noncanonical_DVG_candidate"
            guess_method = "ORF1_alignment_with_noncanonical_structure"
        elif np.isfinite(body_pos):
            best_guess, _ = nearest_sgrna(int(body_pos))
            guess_method = "nearest_body_start_low_confidence"
        else:
            best_guess = "sgRNA_unassigned"
            guess_method = "no_body_alignment"
        method = guess_method
        confidence = "low" if best_guess != "sgRNA_unassigned" else "none"

    largest = max(introns, key=lambda value: value["length"], default=None)
    return {
        "read_id": alignment.query_name,
        "bam_mapq": int(alignment.mapping_quality),
        "bam_query_length_nt": int(alignment.query_length or 0),
        "bam_query_alignment_length_nt": int(alignment.query_alignment_length or 0),
        "reference_start0": int(alignment.reference_start),
        "reference_end0": int(alignment.reference_end),
        "cigar": alignment.cigarstring,
        "n_cigar_N": len(introns),
        "largest_N_length": largest["length"] if largest else np.nan,
        "largest_N_donor0": largest["donor0"] if largest else np.nan,
        "largest_N_acceptor0": largest["acceptor0"] if largest else np.nan,
        "leader_body_candidate_present": bool(chosen is not None),
        "leader_body_donor0": chosen["donor0"] if chosen else np.nan,
        "leader_body_acceptor0": chosen["acceptor0"] if chosen else np.nan,
        "leader_body_intron_length": chosen["length"] if chosen else np.nan,
        "nearest_canonical_sgrna": chosen_label,
        "acceptor_minus_canonical_nt": chosen_delta if chosen is not None else np.nan,
        "junction_left_30nt_aligned": left_flank,
        "junction_right_30nt_aligned": right_flank,
        "leader_body_junction_confirmed": canonical_junction,
        "orf1_unique_aligned_nt": orf1_overlap,
        "grna_orf1_supported": grna_supported,
        "orf1_positive_noncanonical_DVG_candidate": bool(
            orf1_overlap >= args.min_orf1_overlap and len(introns) > 0
        ),
        "leftmost_body_aligned_position0": body_pos,
        "sgrna_primary_assignment": primary,
        "sgrna_all_read_identifiability": (
            primary if primary != "sgRNA_unassigned" else "nested_RNA_unidentifiable"
        ),
        "sgrna_low_confidence_gene_guess": best_guess,
        "sgrna_assignment_method": method,
        "sgrna_assignment_confidence": confidence,
    }


def extract_assignments(
    target: pd.DataFrame,
    bam_dir: Path,
    reference: str,
    args: argparse.Namespace,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for barcode in sorted(target["mapping_barcode"].unique()):
        wanted = set(target.loc[target["mapping_barcode"].eq(barcode), "read_id"])
        bam_path = bam_dir / f"{barcode}.virus.sorted.bam"
        if not bam_path.exists():
            raise FileNotFoundError(bam_path)
        found: set[str] = set()
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for alignment in bam.fetch(reference):
                if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                    continue
                if alignment.query_name not in wanted:
                    continue
                if alignment.query_name in found:
                    raise RuntimeError(f"Duplicate primary alignment: {barcode} {alignment.query_name}")
                row = classify_alignment(alignment, args)
                row["bam_barcode"] = barcode
                row["source_bam"] = str(bam_path)
                records.append(row)
                found.add(alignment.query_name)
        missing = wanted - found
        if missing:
            examples = ",".join(sorted(missing)[:5])
            raise RuntimeError(f"{barcode}: {len(missing)} target reads absent from BAM; examples={examples}")
    assignments = pd.DataFrame(records)
    if assignments["read_id"].duplicated().any():
        raise RuntimeError("read_id is not globally unique in assignment table")
    return assignments
