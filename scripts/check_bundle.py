#!/usr/bin/env python3
"""Check package contents, syntax and citation files."""
import argparse
import ast
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reproduce import ROOT


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--release", action="store_true")
    a = p.parse_args()
    errors, n_python, n_shell = [], 0, 0
    for path in ROOT.rglob("*"):
        rel = path.relative_to(ROOT)
        if any(part in {".git", "runs", "work", "dist", "__pycache__", ".pytest_cache", ".ruff_cache"} for part in rel.parts): continue
        if path.is_symlink(): errors.append(f"Unexpected symlink: {rel}")
        if not path.is_file(): continue
        if path.suffix in {".bam", ".pod5", ".npy", ".parquet", ".blow5", ".fastq", ".pth"} or path.stat().st_size > 10_000_000:
            errors.append(f"Data/large file in source bundle: {rel}")
        if path.suffix == ".py":
            try: ast.parse(path.read_text(), filename=str(path)); n_python += 1
            except SyntaxError as exc: errors.append(f"{rel}: {exc}")
        if path.suffix == ".sh":
            result = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            if result.returncode: errors.append(f"{rel}: {result.stderr}")
            n_shell += 1
    if a.release:
        for name in ("CITATION.cff", "LICENSE"):
            if not (ROOT / name).is_file(): errors.append(f"Release file missing: {name}")
    print(json.dumps({"status": "PASS" if not errors else "FAILED", "python_files": n_python, "shell_files": n_shell, "errors": errors}, indent=2))
    return int(bool(errors))


if __name__ == "__main__":
    raise SystemExit(main())
