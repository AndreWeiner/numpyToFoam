#!/usr/bin/env python3
"""Check imported forces and, when supplied, a native reference."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("reference", type=Path, nargs="?")
    args = parser.parse_args()

    values = np.loadtxt(args.path, comments="#", ndmin=2)
    assert values.shape[1] >= 4, values.shape
    forces = values[:, 1:4]
    assert np.isfinite(forces).all(), forces
    assert np.linalg.norm(forces[-1]) > 0, forces[-1]

    if args.reference:
        reference = np.loadtxt(args.reference, comments="#", ndmin=2)
        np.testing.assert_allclose(values[:, 0], reference[:, 0])
        np.testing.assert_allclose(forces, reference[:, 1:4], rtol=1e-3)

        if values.shape[1] >= 10 and reference.shape[1] >= 10:
            # Pressure is reconstructed from corrected p patches. Turbulent
            # viscosity is regenerated instead of reading stored nut patches.
            np.testing.assert_allclose(
                values[:, 4:7], reference[:, 4:7], rtol=1e-6
            )
            np.testing.assert_allclose(
                values[:, 7:10], reference[:, 7:10], rtol=0.1
            )


if __name__ == "__main__":
    main()
