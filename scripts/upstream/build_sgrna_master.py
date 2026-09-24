#!/usr/bin/env python3
"""Rebuild MHV transcript assignments and read/poly(A)-length metadata."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from transcript_assignment import extract_assignments


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", type=Path, required=True)
    p.add_argument("--anchors", type=Path, required=True)
    p.add_argument("--bam-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    frames = []
    for sample in ("bc04", "bc05", "bc06", "bc07", "bc11"):
        meta = pd.read_parquet(a.features / f"{sample}.metadata.parquet")
        meta = meta.loc[meta.mapping_group.eq("virus"), ["read_id", "mapping_barcode", "warpdemux_polya_estimate_nt", "read_nt_len"]].copy()
        anchors = pd.read_parquet(a.anchors / f"{sample}.dorado_polya_anchors.parquet", columns=["read_id", "dorado_polya_estimate_nt"])
        frames.append(meta.merge(anchors, on="read_id", how="left", validate="one_to_one"))
    target = pd.concat(frames, ignore_index=True)
    if len(target) != 86406 or not target.read_id.is_unique:
        raise RuntimeError("Expected 86,406 unique target MHV reads")
    settings = argparse.Namespace(min_intron=1000, leader_donor_max=100, junction_flank_window=30,
        acceptor_tolerance=30, min_junction_flank_aligned=15, min_mapq=20, min_orf1_overlap=100)
    assigned = extract_assignments(target, a.bam_dir, "NC_048217.1", settings)
    reads = target.merge(assigned, left_on=["read_id", "mapping_barcode"], right_on=["read_id", "bam_barcode"], how="left", validate="one_to_one")
    if reads.sgrna_primary_assignment.isna().any():
        raise RuntimeError("Missing transcript assignment")
    dorado = pd.Series(pd.to_numeric(reads.dorado_polya_estimate_nt, errors="coerce").to_numpy(dtype=float, na_value=np.nan), index=reads.index)
    warp = pd.Series(pd.to_numeric(reads.warpdemux_polya_estimate_nt, errors="coerce").to_numpy(dtype=float, na_value=np.nan), index=reads.index)
    use_dorado = dorado.gt(0) & np.isfinite(dorado)
    reads["polya_length_nt"] = dorado.where(use_dorado, warp)
    reads["polya_length_source"] = np.where(use_dorado, "dorado_pt", "warpdemux_rate_scaled_fallback")
    reads["sgrna_primary_high_confidence"] = reads.sgrna_primary_assignment.ne("sgRNA_unassigned")
    if reads.polya_length_nt.isna().any() or reads.read_nt_len.isna().any():
        raise RuntimeError("Missing length metadata")
    a.output.parent.mkdir(parents=True, exist_ok=True)
    reads.to_parquet(a.output, index=False)
    print(json.dumps({"rows": len(reads), "strict_sgrna": int(reads.sgrna_assignment_method.eq("canonical_leader_body_junction").sum()),
        "dorado_lengths": int(use_dorado.sum()), "settings": vars(settings)}, indent=2))


if __name__ == "__main__":
    main()
