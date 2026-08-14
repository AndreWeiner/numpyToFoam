#!/usr/bin/env python3
"""Validate foamToNumpy function-object output for the cavity test."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


def state_value(path: Path, key: str) -> str:
    match = re.search(rf"^\s*{re.escape(key)}\s+([^;]+);", path.read_text(), re.M)
    if not match:
        raise AssertionError(f"Missing {key!r} in {path}")
    return match.group(1).strip()


def validate(root: Path, expected_times: list[float]) -> None:
    segments = sorted(path for path in root.iterdir() if path.is_dir())
    assert len(segments) == 1, segments
    segment = segments[0]

    batches = sorted(segment.glob("batch_*"))
    expected_counts = [2] * (len(expected_times) // 2)
    if len(expected_times) % 2:
        expected_counts.append(1)

    assert len(batches) == len(expected_counts), batches
    loaded_times: list[float] = []

    for batch, expected_count in zip(batches, expected_counts):
        state = batch / "state"
        count = int(state_value(state, "count"))
        assert count == expected_count
        assert state_value(state, "sealed") == "true"

        times = np.load(batch / "times.npy", allow_pickle=False)
        assert times.shape == (count,)
        assert times.flags.f_contiguous
        loaded_times.extend(times.tolist())

        for field_name, components in (("p", 1), ("U", 3)):
            files = sorted(batch.glob(f"{field_name}_proc_*.npy"))
            assert files, (batch, field_name)

            for path in files:
                array = np.load(path, allow_pickle=False)
                expected_shape = (
                    (array.shape[0], count)
                    if components == 1
                    else (array.shape[0], components, count)
                )
                assert array.shape == expected_shape
                assert array.dtype == np.float64
                assert array.flags.f_contiguous

    np.testing.assert_allclose(loaded_times, expected_times)

    geometry = sorted(segment.glob("geometry_*"))
    assert len(geometry) == 1, geometry
    for pattern in ("cellCentres_proc_*.npy", "cellVolumes_proc_*.npy"):
        files = sorted(geometry[0].glob(pattern))
        assert files
        for path in files:
            array = np.load(path, allow_pickle=False)
            assert array.shape[-1] == 1
            assert array.flags.f_contiguous


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("times", nargs="+", type=float)
    args = parser.parse_args()
    validate(args.root, args.times)


if __name__ == "__main__":
    main()
