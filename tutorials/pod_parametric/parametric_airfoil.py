#!/usr/bin/env python3
"""Small parametric POD interpolation example for the airfoil case."""

import csv
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


CASE = Path("of_airfoil")
PODI_CASE = Path("of_airfoil_podi")
PREDICTIONS = Path("prediction_data")
RESULTS = Path("results")
BASE_SPEED = np.linalg.norm([25.75, 3.62])
INITIAL_SPEED = 21.777687887
N_MODES = 5

# The first six samples train a quadratic response surface; the last two test it.
SAMPLES = [
    (-10.0, 0.75),
    (-10.0, 1.25),
    (10.0, 0.75),
    (10.0, 1.25),
    (0.0, 1.00),
    (-5.0, 1.00),
    (5.0, 1.00),
    (0.0, 0.875),
]


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


def run(*command, cwd=None):
    print("+", *command)
    subprocess.run(command, cwd=cwd, check=True)


def configure(alpha_deg, speed):
    alpha = math.radians(alpha_deg)
    drag = np.array([math.cos(alpha), math.sin(alpha), 0.0])
    lift = np.array([-math.sin(alpha), math.cos(alpha), 0.0])
    velocity = speed * drag

    U_path = CASE / "0.orig" / "U"
    text = U_path.read_text()
    text = re.sub(
        r"internalField\s+uniform\s+\([^;]+\);",
        f"internalField   uniform ({velocity[0]:.12g} {velocity[1]:.12g} 0);",
        text,
        count=1,
    )
    U_path.write_text(text)

    control_path = CASE / "system" / "controlDict"
    text = control_path.read_text()
    replacements = {
        r"magUInf\s+[-+0-9.eE]+;": f"magUInf         {speed:.12g};",
        r"liftDir\s+\([^;]+\);":
            f"liftDir         ({lift[0]:.12g} {lift[1]:.12g} 0);",
        r"dragDir\s+\([^;]+\);":
            f"dragDir         ({drag[0]:.12g} {drag[1]:.12g} 0);",
    }
    for pattern, replacement in replacements.items():
        text = re.sub(pattern, replacement, text, count=1)
    control_path.write_text(text)


def load_export():
    segment = CASE / "postProcessing" / "numpyExport" / "0"
    batch = segment / "batch_000000"
    p_parts = [np.load(path)[:, 0] for path in sorted(batch.glob("p_proc_*.npy"))]
    U_parts = [np.load(path)[:, :, 0] for path in sorted(batch.glob("U_proc_*.npy"))]
    nu_tilda_parts = [
        np.load(path)[:, 0]
        for path in sorted(batch.glob("nuTilda_proc_*.npy"))
    ]
    volumes = np.concatenate([
        np.load(path)[:, 0]
        for path in sorted((segment / "geometry_000000").glob(
            "cellVolumes_proc_*.npy"
        ))
    ])
    counts = [len(part) for part in p_parts]
    p = np.concatenate(p_parts)
    U = np.concatenate(U_parts)
    nu_tilda = np.concatenate(nu_tilda_parts)
    return np.concatenate((p, U[:, 0], U[:, 1], nu_tilda)), counts, volumes


def force_coefficients():
    path = sorted((CASE / "postProcessing" / "forceCoeffs").glob("**/*.dat"))[-1]
    row = np.loadtxt(path, comments="#", ndmin=2)[-1]
    return row[[1, 3, 5]]


def basis(parameters):
    alpha = parameters[:, 0] / 10.0
    speed = (parameters[:, 1] / BASE_SPEED - 1.0) / 0.25
    return np.column_stack((
        np.ones(len(parameters)),
        alpha,
        speed,
        alpha * speed,
        alpha**2,
        speed**2,
    ))


def write_predictions(state, counts):
    shutil.rmtree(PREDICTIONS, ignore_errors=True)
    batch = PREDICTIONS / "0" / "batch_000000"
    batch.mkdir(parents=True)
    times = np.arange(1, state.shape[1] + 1, dtype=float)
    save_fortran(batch / "times.npy", times)

    n_cells = sum(counts)
    p = state[:n_cells]
    U = np.zeros((n_cells, 3, state.shape[1]), order="F")
    U[:, 0] = state[n_cells:2 * n_cells]
    U[:, 1] = state[2 * n_cells:3 * n_cells]
    nu_tilda = state[3 * n_cells:]

    start = 0
    for proci, count in enumerate(counts):
        stop = start + count
        save_fortran(batch / f"p_proc_{proci}.npy", p[start:stop])
        save_fortran(batch / f"U_proc_{proci}.npy", U[start:stop])
        save_fortran(
            batch / f"nuTilda_proc_{proci}.npy",
            nu_tilda[start:stop],
        )
        start = stop

    (PREDICTIONS / "0" / "segmentInfo").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object segmentInfo; }\n"
    )
    (batch / "state").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object state; }\n"
        f"count {len(times)};\nmeshRevision 0;\nsealed true;\n"
    )


def copy_case():
    def ignore(_, names):
        return {name for name in names
                if name in {"500", "postProcessing"} or name.startswith("log.")}

    shutil.rmtree(PODI_CASE, ignore_errors=True)
    shutil.copytree(CASE, PODI_CASE, ignore=ignore)
    shutil.copytree(PREDICTIONS, PODI_CASE / PREDICTIONS)


def imported_coefficients(parameters):
    output = PODI_CASE / "postProcessing" / "importedForces"
    force_path = sorted(
        output.glob("**/force.dat")
    )[-1]
    moment_path = sorted(output.glob("**/moment.dat"))[-1]
    forces = np.loadtxt(force_path, comments="#", ndmin=2)[:, 1:4]
    moments = np.loadtxt(moment_path, comments="#", ndmin=2)[:, 1:4]
    coefficients = []

    for force, moment, (alpha_deg, speed) in zip(
        forces, moments, parameters
    ):
        alpha = math.radians(alpha_deg)
        drag = np.array([math.cos(alpha), math.sin(alpha), 0.0])
        lift = np.array([-math.sin(alpha), math.cos(alpha), 0.0])
        q = 0.5 * speed**2
        coefficients.append([
            np.dot(force, drag) / q,
            np.dot(force, lift) / q,
            moment[2] / q,
        ])

    return np.asarray(coefficients)


shutil.rmtree(RESULTS, ignore_errors=True)
RESULTS.mkdir()
parameters = np.array([(alpha, factor * BASE_SPEED) for alpha, factor in SAMPLES])
snapshots = []
reference_coefficients = []

for casei, (alpha, speed) in enumerate(parameters):
    print(f"\nCase {casei}: alpha={alpha:g}, speed={speed:g}")
    configure(alpha, speed)
    run("./Allclean", cwd=CASE)
    shutil.rmtree(CASE / "postProcessing", ignore_errors=True)
    run("./Allrun", cwd=CASE)
    snapshot, counts, volumes = load_export()
    snapshots.append(snapshot)
    reference_coefficients.append(force_coefficients())

state = np.column_stack(snapshots)
n_cells = sum(counts)
scale = np.concatenate((
    np.full(n_cells, BASE_SPEED**2),
    np.full(2 * n_cells, BASE_SPEED),
    np.full(n_cells, 4e-5),
))[:, None]
train = np.arange(6)
mean = (state / scale)[:, train].mean(axis=1, keepdims=True)
weights = np.sqrt(np.tile(volumes, 4))[:, None]
modes, singular_values, _ = np.linalg.svd(
    weights * ((state / scale)[:, train] - mean),
    full_matrices=False,
)
modes = modes[:, :N_MODES]
coefficients = modes.T @ (weights * ((state / scale)[:, train] - mean))
response = np.linalg.lstsq(basis(parameters[train]), coefficients.T, rcond=None)[0]
prediction = scale * (mean + modes @ (basis(parameters) @ response).T / weights)

write_predictions(prediction, counts)
copy_case()
configure(-6.5, INITIAL_SPEED)
run(
    "mpirun",
    "-np",
    str(len(counts)),
    "numpyPostProcess",
    "-parallel",
    cwd=PODI_CASE,
)

predicted_coefficients = imported_coefficients(parameters)
with (RESULTS / "comparison.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow((
        "case", "training", "alpha", "speed",
        "Cd", "Cd_PODI", "Cl", "Cl_PODI", "Cm", "Cm_PODI",
    ))
    for casei, ((alpha, speed), reference, predicted) in enumerate(zip(
        parameters,
        reference_coefficients,
        predicted_coefficients,
    )):
        writer.writerow((
            casei, casei in train, alpha, speed,
            reference[0], predicted[0],
            reference[1], predicted[1],
            reference[2], predicted[2],
        ))

np.save(RESULTS / "singular_values.npy", singular_values)
print(f"Wrote {RESULTS / 'comparison.csv'}")
