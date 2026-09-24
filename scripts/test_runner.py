"""Tests for runner planning, output paths and stage execution."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import reproduce as r


class RunnerSafetyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def plan(self, out):
        return {"target": "publication_terminal_nona_results_v1", "output": str(out), "project_root": str(self.root),
                "rebuild_inputs": False, "required_inputs": [], "steps": []}

    def test_existing_output_is_refused_and_untouched(self):
        out = self.root / "existing"
        out.mkdir()
        marker = out / "keep.txt"
        marker.write_text("original")
        self.assertTrue(any("Output already exists" in e for e in r.preflight(self.plan(out))))
        self.assertEqual(marker.read_text(), "original")

    def test_dangling_link_is_refused(self):
        out = self.root / "link"
        out.symlink_to(self.root / "missing")
        self.assertTrue(any("Output already exists" in e for e in r.preflight(self.plan(out))))
        self.assertFalse((self.root / "missing").exists())

    def test_missing_input_does_not_create_output(self):
        plan = self.plan(self.root / "out")
        plan["required_inputs"] = [str(self.root / "missing.pod5")]
        self.assertTrue(any("Missing/empty" in e for e in r.preflight(plan)))
        self.assertFalse((self.root / "out").exists())

    def test_plan_is_read_only_and_only_uses_target_method(self):
        out = self.root / "out"
        with contextlib.redirect_stdout(io.StringIO()) as stream:
            rc = r.main(["--project-root", str(self.root / "no_data"), "--output", str(out), "--plan"])
        plan = json.loads(stream.getvalue())
        self.assertEqual(rc, 0)
        self.assertEqual(plan["target"], "publication_terminal_nona_results_v1")
        self.assertEqual(plan["steps"][0]["name"], "01_mapping_wdx_metadata")
        self.assertEqual(plan["analysis"], "barcode_figure")
        self.assertFalse(out.exists())

    def test_failed_process_stops_following_steps(self):
        out = self.root / "out"
        plan = self.plan(out)
        plan["steps"] = [{"name": "first", "argv": ["fake-first"]}, {"name": "second", "argv": ["fake-second"]}]
        with patch.object(r.subprocess, "run") as proc, contextlib.redirect_stdout(io.StringIO()):
            proc.return_value.returncode = 7
            with self.assertRaises(RuntimeError): r.execute(plan)
            self.assertEqual(proc.call_count, 1)
        report = json.loads((out / "run_manifest.json").read_text())
        self.assertEqual(report["status"], "FAILED")
        self.assertEqual(report["steps"][0]["returncode"], 7)

    def test_rebuild_does_not_require_precomputed_parquet_or_reference_results(self):
        with contextlib.redirect_stdout(io.StringIO()) as stream:
            rc = r.main(["--project-root", str(self.root / "data"), "--output", str(self.root / "out"),
                         "--rebuild-inputs", "--plan"])
        plan = json.loads(stream.getvalue())
        self.assertEqual(rc, 0)
        self.assertEqual(len(plan["steps"]), 6)
        self.assertFalse(any(name.endswith(".parquet") or "/results/" in name for name in plan["required_inputs"]))
        self.assertFalse(any(arg.startswith("--reference-") for step in plan["steps"] for arg in step["argv"]))
        self.assertFalse(any(name.endswith(".bai") for name in plan["required_inputs"]))
        self.assertFalse(any("sgrna" in arg for step in plan["steps"] for arg in step["argv"]))

    def test_advanced_workflow_includes_transcript_analysis(self):
        with contextlib.redirect_stdout(io.StringIO()) as stream:
            r.main(["--project-root", str(self.root / "data"), "--output", str(self.root / "out"),
                    "--advanced", "--plan"])
        plan = json.loads(stream.getvalue())
        self.assertEqual(plan["analysis"], "advanced")
        self.assertEqual(len(plan["steps"]), 11)
        self.assertTrue(any("build_sgrna_master.py" in arg for step in plan["steps"] for arg in step["argv"]))

    def test_figure_cached_inputs_need_signal_and_dorado_metadata(self):
        with contextlib.redirect_stdout(io.StringIO()) as stream:
            r.main(["--project-root", str(self.root / "data"), "--output", str(self.root / "out"),
                    "--inputs", str(self.root / "cache"), "--plan"])
        plan = json.loads(stream.getvalue())
        self.assertFalse(plan["rebuild_inputs"])
        self.assertEqual(len(plan["required_inputs"]), 15)
        self.assertFalse(any("sgrna" in name for name in plan["required_inputs"]))


if __name__ == "__main__":
    unittest.main()
