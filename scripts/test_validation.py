"""Regression checks for direct numerical result validation."""
import unittest

import pandas as pd

from validate_figure import check_calls, check_table, compare_reads


class ResultValidationTests(unittest.TestCase):
    def test_direct_comparison_detects_a_changed_score(self):
        expected = pd.DataFrame({"read_id": ["a", "b"], "score": [1.0, 2.0]})
        actual = expected.copy()
        actual.loc[1, "score"] = 2.000001
        with self.assertRaises(AssertionError):
            compare_reads(actual, expected)

    def test_direct_comparison_accepts_reordered_reads(self):
        expected = pd.DataFrame({"read_id": ["a", "b"], "score": [1.0, 2.0]})
        compare_reads(expected.iloc[::-1], expected)

    def test_duplicate_read_ids_are_rejected(self):
        frame = pd.DataFrame({"read_id": ["a", "a"], "score": [1.0, 2.0]})
        with self.assertRaises(ValueError):
            check_table(frame, 2, ["read_id", "score"])

    def test_call_validation_checks_cutoff_and_missing_calls(self):
        frame = pd.DataFrame({
            "eligible_terminal_nona_direct": [1, 0],
            "terminal_nona_signed_difference_pa": [-5.0, float("nan")],
            "terminal_nona_abs_difference_pa": [5.0, float("nan")],
            "terminal_junction_nonA_like": pd.Series([True, pd.NA], dtype="boolean"),
        })
        check_calls(frame, 5.0)
        frame.loc[0, "terminal_junction_nonA_like"] = False
        with self.assertRaises(ValueError):
            check_calls(frame, 5.0)
        frame.loc[0, "terminal_junction_nonA_like"] = True
        frame.loc[1, "terminal_junction_nonA_like"] = False
        with self.assertRaises(ValueError):
            check_calls(frame, 5.0)


if __name__ == "__main__":
    unittest.main()
