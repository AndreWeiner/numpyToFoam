#!/usr/bin/env python3
"""Compute and write a weighted POD reconstruction and its first five modes."""

from pathlib import Path
import shutil

import numpy as np


SOURCE = Path("of_cavity/postProcessing/numpyExport/0")
RANK = 15
N_MODES = 5


def save_fortran(path, values):
    values = np.asfortranarray(values)
    header = {
        "descr": np.lib.format.dtype_to_descr(values.dtype),
        "fortran_order": True,
        "shape": values.shape,
    }
    with path.open("wb") as stream:
        np.lib.format.write_array_header_1_0(stream, header)
        stream.write(values.tobytes(order="F"))


def load_parts(field):
    batches = sorted(SOURCE.glob("batch_*"))
    names = sorted(path.name for path in batches[0].glob(f"{field}_proc_*.npy"))
    return [np.concatenate([np.load(batch / name) for batch in batches], axis=-1)
            for name in names]


def write_dataset(path, times, p, U, counts):
    shutil.rmtree(path, ignore_errors=True)
    batch = path / "0" / "batch_000000"
    batch.mkdir(parents=True)

    save_fortran(batch / "times.npy", times)
    start = 0
    for proci, count in enumerate(counts):
        stop = start + count
        save_fortran(batch / f"p_proc_{proci}.npy", p[start:stop])
        save_fortran(batch / f"U_proc_{proci}.npy", U[start:stop])
        start = stop

    (path / "0" / "segmentInfo").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object segmentInfo; }\n"
    )
    (batch / "state").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object state; }\n"
        f"count {len(times)};\nmeshRevision 0;\nsealed true;\n"
    )


p_parts = load_parts("p")
U_parts = load_parts("U")
counts = [part.shape[0] for part in p_parts]
p = np.concatenate(p_parts)
U = np.concatenate(U_parts)
volumes = np.concatenate([
    np.load(path)[:, 0]
    for path in sorted((SOURCE / "geometry_000000").glob("cellVolumes_proc_*.npy"))
])
times = np.concatenate([
    np.load(path)
    for path in sorted(SOURCE.glob("batch_*/times.npy"))
])

n_cells = p.shape[0]
state = np.vstack((p, U[:, 0], U[:, 1]))
weights = np.sqrt(np.tile(volumes, 3))[:, None]
modes, singular_values, coefficients = np.linalg.svd(
    weights * state,
    full_matrices=False,
)

rank = min(RANK, len(singular_values))
reconstruction = (
    modes[:, :rank]
    @ np.diag(singular_values[:rank])
    @ coefficients[:rank]
) / weights
selected_modes = modes[:, :N_MODES] / weights


def split_state(values):
    velocity = np.zeros((n_cells, 3, values.shape[1]), order="F")
    velocity[:, 0] = values[n_cells:2 * n_cells]
    velocity[:, 1] = values[2 * n_cells:]
    return values[:n_cells], velocity


p_reconstructed, U_reconstructed = split_state(reconstruction)
p_modes, U_modes = split_state(selected_modes)

write_dataset(
    Path("reconstruction_data"),
    times,
    p_reconstructed,
    U_reconstructed,
    counts,
)
write_dataset(
    Path("mode_data"),
    np.arange(1, N_MODES + 1, dtype=float),
    p_modes,
    U_modes,
    counts,
)

print(f"Retained {rank} modes for the reconstruction")
