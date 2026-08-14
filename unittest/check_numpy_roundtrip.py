#!/usr/bin/env python3
"""Check byte-identical NumPy export/import/re-export snapshots."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
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
            assert target.is_file(), target

            source_array = np.load(source, allow_pickle=False)
            target_array = np.load(target, allow_pickle=False)
            assert source_array.flags.f_contiguous, source
            assert target_array.flags.f_contiguous, target
            np.testing.assert_array_equal(target_array, source_array)
            assert digest(target) == digest(source), (source, target)


if __name__ == "__main__":
    main()
