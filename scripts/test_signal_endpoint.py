"""Scientific boundary checks for window geometry, QC and the fixed cutoff."""
from argparse import Namespace
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import call_terminal_nona_direct_cutoff as caller


def arguments(**kwargs):
    return Namespace(anchor_mode="warp_raw", boundary_method="cnn", cutoff_pa=5.0,
                     right_phase_nt=0.5, fixed_rna_dwell_samples=4000 / 130,
                     reference_min_polya_nt=100.0, reference_min_qscore=12.0, **kwargs)


class SignalEndpointTests(unittest.TestCase):
    def test_raw_window_positions_and_independent_callability_gates(self):
        ids = ["complete", "short_tail", "fallback", "left_edge", "nonfinite", "missing"]
        frame = pd.DataFrame({
            "read_id": ids, "polya_start": [20, 20, 20, 15, 20, 20],
            "polya_end": [120, 96, 120, 120, 120, 120],
            "mapping_group": "virus", "mapping_barcode": "bc05",
            "barcode_mapping_consistent": 1, "boundary_method": ["cnn", "cnn", "fallback", "cnn", "cnn", "cnn"],
            "dorado_polya_start": np.nan, "dorado_pa_anchor": np.nan,
        })
        records = []
        for name in ids[:-1]:
            signal = np.arange(120, dtype=np.float32)
            if name == "nonfinite":
                signal[70] = np.nan
            records.append(SimpleNamespace(read_id=name, signal_pa=signal))
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / "bc05.pod5").touch()
            with patch("pod5.Reader") as reader:
                reader.return_value.__enter__.return_value.reads.return_value = records
                scored = caller.prepare_barcode_raw("bc05", frame, arguments(pod5_dir=Path(folder)))
        row = scored.iloc[0]
        self.assertEqual(row["terminal_nona_terminal_left_array_index"], 4)
        self.assertEqual(row["terminal_nona_model_right_array_index"], 35)
        self.assertEqual(row["terminal_window_mean_centered_pa"], 19)
        self.assertEqual(row["stableA_window_mean_centered_pa"], 81)
        self.assertEqual(row["terminal_minus_stableA_pa"], -62)
        self.assertEqual(scored["eligible_terminal_nona_direct"].tolist(), [1, 0, 0, 0, 0, 0])
        self.assertEqual(scored["terminal_nona_direct_qc_status"].tolist(), [
            "eligible", "warp_polya_too_short_for_stableA_window", "nonprimary_boundary_method",
            "terminal_or_stableA_window_incomplete", "terminal_or_stableA_window_incomplete", "pod5_signal_not_found",
        ])

    def test_two_sided_inclusive_cutoff_and_missing_calls(self):
        args = arguments()
        host = pd.DataFrame({
            "read_id": [f"host-{i}" for i in range(1002)], "mapping_group": "host",
            "eligible_terminal_nona_direct": 1, "actual_polya_len_samples": 4000,
            "dorado_mean_qscore": 12.0, "terminal_minus_stableA_pa": np.float32(0),
        })
        host.loc[1000, "dorado_mean_qscore"] = 11.99
        host.loc[1001, "actual_polya_len_samples"] = 3000
        just_below = np.nextafter(np.float32(5), np.float32(0))
        virus = pd.DataFrame({
            "read_id": [f"virus-{i}" for i in range(5)], "mapping_group": "virus",
            "eligible_terminal_nona_direct": [1, 1, 1, 1, 0],
            "actual_polya_len_samples": 4000, "dorado_mean_qscore": 12.0,
            "terminal_minus_stableA_pa": np.array([-5, 5, -just_below, just_below, 10], dtype=np.float32),
        })
        frame = pd.concat([host, virus], ignore_index=True)
        with patch.object(caller, "BARCODES", ("bc05",)), patch.object(caller, "prepare_barcode", return_value=frame):
            scored, membership, reference, _ = caller.score_all_reads(args)
        self.assertEqual(reference, 0)
        self.assertEqual(int(membership.sum()), 1000)
        observed = scored.loc[scored["mapping_group"].eq("virus")]
        self.assertEqual(observed["terminal_junction_nonA_like"].iloc[:4].tolist(), [True, True, False, False])
        self.assertTrue(pd.isna(observed.iloc[4]["terminal_junction_nonA_like"]))
        self.assertTrue(pd.isna(observed.iloc[4]["terminal_nona_abs_difference_pa"]))
        self.assertEqual(observed.iloc[4]["terminal_nona_direct_class"], "uncallable")


if __name__ == "__main__":
    unittest.main()
