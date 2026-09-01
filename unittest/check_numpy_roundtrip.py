#!/usr/bin/env python3
"""Check byte-identical NumPy export/import/re-export snapshots."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_arrays(source: Path, target: Path) -> np.ndarray:
    assert target.is_file(), target

    source_array = np.load(source, allow_pickle=False)
    target_array = np.load(target, allow_pickle=False)
    assert source_array.flags.f_contiguous, source
    assert target_array.flags.f_contiguous, target
    np.testing.assert_array_equal(target_array, source_array)
    assert digest(target) == digest(source), (source, target)
    return source_array


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", action="store_true")
    parser.add_argument("source", type=Path)
    parser.add_argument("reexport", type=Path)
    parser.add_argument("fields", nargs="+")
    args = parser.parse_args()

    for field in [*args.fields, "times"]:
        source_files = sorted(args.source.glob(f"batch_*/{field}_proc_*.npy"))
        if field == "times":
            source_files = sorted(args.source.glob("batch_*/times.npy"))

        assert source_files, (args.source, field)

        for source in source_files:
            relative = source.relative_to(args.source)
            target = args.reexport / relative
            compare_arrays(source, target)

    if args.geometry:
        geometry_files = sorted(args.source.glob("geometry_*/*.npy"))
        assert geometry_files, args.source
        assert any(path.name.startswith("faceAreas_proc_") for path in geometry_files)

        for source in geometry_files:
            relative = source.relative_to(args.source)
            source_array = compare_arrays(source, args.reexport / relative)

            if source.name.startswith("faceAreas_proc_"):
                assert source_array.ndim == 2, source
                assert source_array.shape[-1] == 1, source
                assert np.all(source_array > 0), source


if __name__ == "__main__":
    main()
