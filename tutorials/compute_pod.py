#!/usr/bin/env python3
"""Compute POD data from pressure and in-plane velocity snapshots."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import numpy as np

PRESSURE_FIELD = "p"
VELOCITY_FIELD = "U"
PRESSURE_MODE_FIELD = "p_modes"
VELOCITY_MODE_FIELD = "U_modes"
VOLUME_FIELD = "cellVolumes"


def processor_files(data_dir: Path, field: str) -> list[tuple[int, Path]]:
    field_dir = data_dir / field
    pattern = re.compile(rf"^{re.escape(field)}_proc_(\d+)\.npy$")
    matches: list[tuple[int, Path]] = []

    for path in field_dir.glob(f"{field}_proc_*.npy"):
        match = pattern.match(path.name)
        if match:
            matches.append((int(match.group(1)), path))

    if not matches:
        raise FileNotFoundError(f"No processor files found in {field_dir}")

    matches.sort(key=lambda item: item[0])
    return matches


def load_scalar_parts(data_dir: Path, field: str) -> tuple[list[int], list[np.ndarray]]:
    proc_ids: list[int] = []
    arrays: list[np.ndarray] = []

    for proc, path in processor_files(data_dir, field):
        array = np.load(path)
        if array.ndim != 2:
            raise ValueError(
                f"Expected scalar field shape (nCells, nTimes), got "
                f"{array.shape} in {path}"
            )
        proc_ids.append(proc)
        arrays.append(np.asarray(array, dtype=np.float64))

    return proc_ids, arrays


def load_vector_parts(data_dir: Path, field: str) -> tuple[list[int], list[np.ndarray]]:
    proc_ids: list[int] = []
    arrays: list[np.ndarray] = []

    for proc, path in processor_files(data_dir, field):
        array = np.load(path)
        if array.ndim != 3 or array.shape[1] < 2:
            raise ValueError(
                f"Expected vector field shape (nCells, 3, nTimes), got "
                f"{array.shape} in {path}"
            )
        proc_ids.append(proc)
        arrays.append(np.asarray(array, dtype=np.float64))

    return proc_ids, arrays


def load_volume_parts(data_dir: Path) -> tuple[list[int], list[np.ndarray]]:
    proc_ids, arrays = load_scalar_parts(data_dir, VOLUME_FIELD)
    volumes = [array.reshape(-1) for array in arrays]

    for array in volumes:
        if np.any(array <= 0.0):
            raise ValueError("Cell volumes must be positive")

    return proc_ids, volumes


def reset_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)


def write_scalar_parts(
    data: np.ndarray,
    proc_ids: list[int],
    counts: list[int],
    output_dir: Path,
    field: str,
) -> None:
    field_dir = output_dir / field
    field_dir.mkdir(parents=True, exist_ok=True)

    start = 0
    for proc, count in zip(proc_ids, counts):
        stop = start + count
        np.save(field_dir / f"{field}_proc_{proc}.npy", np.asfortranarray(data[start:stop]))
        start = stop


def write_vector_parts(
    x_data: np.ndarray,
    y_data: np.ndarray,
    proc_ids: list[int],
    counts: list[int],
    output_dir: Path,
    field: str,
) -> None:
    field_dir = output_dir / field
    field_dir.mkdir(parents=True, exist_ok=True)

    start = 0
    for proc, count in zip(proc_ids, counts):
        stop = start + count
        part = np.zeros((count, 3, x_data.shape[1]), dtype=np.float64, order="F")
        part[:, 0, :] = x_data[start:stop]
        part[:, 1, :] = y_data[start:stop]
        np.save(field_dir / f"{field}_proc_{proc}.npy", part)
        start = stop


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute POD reconstruction and modes from p, Ux, and Uy "
            "foamToNumpy snapshots."
        )
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--reconstruction-output", type=Path, required=True)
    parser.add_argument("--modes-output", type=Path, required=True)
    parser.add_argument("--rank", type=int, default=15)
    parser.add_argument("--n-modes", type=int, default=5)
    args = parser.parse_args()

    p_proc_ids, p_parts = load_scalar_parts(args.input, PRESSURE_FIELD)
    u_proc_ids, u_parts = load_vector_parts(args.input, VELOCITY_FIELD)
    vol_proc_ids, volume_parts = load_volume_parts(args.input)

    if p_proc_ids != u_proc_ids or p_proc_ids != vol_proc_ids:
        raise ValueError(
            "Processor ids differ between exported fields: "
            f"{PRESSURE_FIELD}={p_proc_ids}, "
            f"{VELOCITY_FIELD}={u_proc_ids}, "
            f"{VOLUME_FIELD}={vol_proc_ids}"
        )

    counts = [part.shape[0] for part in p_parts]
    if counts != [part.shape[0] for part in u_parts]:
        raise ValueError("Pressure and velocity cell counts differ")
    if counts != [part.size for part in volume_parts]:
        raise ValueError("Field and cell-volume counts differ")

    p_snap = np.concatenate(p_parts, axis=0)
    u_snap = np.concatenate(u_parts, axis=0)
    volumes = np.concatenate(volume_parts)

    if p_snap.shape[1] != u_snap.shape[2]:
        raise ValueError(
            f"Pressure has {p_snap.shape[1]} times, velocity has {u_snap.shape[2]}"
        )

    n_cells = p_snap.shape[0]
    state = np.vstack((p_snap, u_snap[:, 0, :], u_snap[:, 1, :]))
    weights = np.sqrt(np.tile(volumes, 3))[:, np.newaxis]
    weighted_state = weights * state

    modes, singular_values, temporal_coefficients = np.linalg.svd(
        weighted_state, full_matrices=False
    )

    rank = min(args.rank, singular_values.size)
    n_modes = min(args.n_modes, modes.shape[1])

    weighted_reconstruction = (
        modes[:, :rank]
        @ np.diag(singular_values[:rank])
        @ temporal_coefficients[:rank, :]
    )
    reconstruction = weighted_reconstruction / weights
    selected_modes = modes[:, :n_modes] / weights

    reset_dir(args.reconstruction_output)
    reset_dir(args.modes_output)

    p_reconstruction = reconstruction[:n_cells]
    ux_reconstruction = reconstruction[n_cells : 2 * n_cells]
    uy_reconstruction = reconstruction[2 * n_cells :]

    p_modes = selected_modes[:n_cells]
    ux_modes = selected_modes[n_cells : 2 * n_cells]
    uy_modes = selected_modes[2 * n_cells :]

    write_scalar_parts(
        p_reconstruction,
        p_proc_ids,
        counts,
        args.reconstruction_output,
        PRESSURE_FIELD,
    )
    write_vector_parts(
        ux_reconstruction,
        uy_reconstruction,
        p_proc_ids,
        counts,
        args.reconstruction_output,
        VELOCITY_FIELD,
    )
    write_scalar_parts(
        p_modes,
        p_proc_ids,
        counts,
        args.modes_output,
        PRESSURE_MODE_FIELD,
    )
    write_vector_parts(
        ux_modes,
        uy_modes,
        p_proc_ids,
        counts,
        args.modes_output,
        VELOCITY_MODE_FIELD,
    )

    print(
        f"Wrote rank-{rank} reconstruction for weighted state shape "
        f"{weighted_state.shape} to {args.reconstruction_output}"
    )
    print(
        f"Wrote {n_modes} POD modes for pressure and velocity "
        f"to {args.modes_output}"
    )


if __name__ == "__main__":
    main()
