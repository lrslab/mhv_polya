"""Verify that barcode count labels follow the plotted source data."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib.figure
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import build_figure6A_overall_mut_vs_wt as plot


class FigureLabelTests(unittest.TestCase):
    def test_count_annotation_changes_with_summary(self):
        root = Path(__file__).resolve().parents[1]
        summary = pd.read_csv(root / "expected/figure6A/figure6A_overall_virus_barcode_summary.tsv", sep="\t")
        row = summary["mapping_barcode"].eq("bc05")
        summary.loc[row, "n_terminal_junction_nonA_like"] = 1202
        summary.loc[row, "terminal_junction_nonA_like_fraction"] = 1202 / 9370
        summary.loc[row, "A_junction_like_fraction"] = 1 - 1202 / 9370
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(matplotlib.figure.Figure, "savefig", autospec=True) as save:
                plot.draw_barcode_figure(summary, Path(folder), dpi=100)
                figure = save.call_args.args[0]
                labels = [text.get_text() for axis in figure.axes for text in axis.texts]
        self.assertTrue(any("WT / bc05: 1,202 / 9,370 (12.83%)" in label for label in labels))
        self.assertFalse(any("WT / bc05: 1,201" in label for label in labels))


if __name__ == "__main__":
    unittest.main()
