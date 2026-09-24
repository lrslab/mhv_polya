#!/usr/bin/env python3
"""Package the code, documentation and expected result tables in a ZIP."""
import argparse
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE_DIRS = {".git", "runs", "work", "dist", "__pycache__", ".pytest_cache", ".ruff_cache"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    output = a.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    subprocess.run([sys.executable, str(ROOT / "scripts/check_bundle.py")], check=True)
    paths = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if not path.is_file() or path.is_symlink() or any(part in EXCLUDE_DIRS for part in relative.parts):
            continue
        if path == output or path.name == ".DS_Store":
            continue
        if path.name.startswith("remote_") and path.name.endswith(".tar.gz"):
            continue
        paths.append(path)
    paths.sort(key=lambda path: str(path.relative_to(ROOT)))
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(paths):
            info = zipfile.ZipInfo("mhv_polya/" + str(path.relative_to(ROOT)), date_time=(2026, 9, 23, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o100755 if path.suffix == ".sh" else 0o100644) << 16
            archive.writestr(info, path.read_bytes())
    with zipfile.ZipFile(output) as archive:
        bad = archive.testzip()
        if bad:
            raise RuntimeError(f"Invalid archive entry: {bad}")
    print(f"Created {output}: {len(paths)} files; {output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
