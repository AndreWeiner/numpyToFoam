#!/usr/bin/env python3
"""Parametric airfoil POD with interpolation."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

CASE_DIR = Path("of_airfoil")
RESULTS_DIR = Path("results")
PREDICTION_DATA = Path("prediction_data")
PODI_CASE = Path("of_airfoil_podi")

BASE_U = np.array([25.75, 3.62], dtype=float)
BASE_SPEED = float(np.linalg.norm(BASE_U))
REF_SPEED = BASE_SPEED
REF_PRESSURE = REF_SPEED**2
ALPHA_RANGE = (-10.0, 10.0)
SPEED_RANGE = (0.75 * BASE_SPEED, 1.25 * BASE_SPEED)
MIN_CASES = 7


def run(command: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


def sobol_2d(n_points: int, bits: int = 30) -> np.ndarray:
    """Small two-dimensional Sobol sequence, enough for this tutorial."""
    dirs_x = [1 << (bits - j) for j in range(1, bits + 1)]
    dirs_y = [1 << (bits - 1)]
    for _ in range(1, bits):
        previous = dirs_y[-1]
        dirs_y.append(previous ^ (previous >> 1))

    points = np.empty((n_points, 2), dtype=float)
    scale = float(1 << bits)
    for i in range(n_points):
        index = i + 1
        x = 0
        y = 0
        bit = 0
        while index:
            if index & 1:
                x ^= dirs_x[bit]
                y ^= dirs_y[bit]
            index >>= 1
            bit += 1
        points[i] = (x / scale, y / scale)
    return points


def latinized_sobol(n_points: int) -> np.ndarray:
    sobol = sobol_2d(n_points)
    design = np.empty_like(sobol)
    for dim in range(2):
        ranks = np.argsort(np.argsort(sobol[:, dim]))
        design[:, dim] = (ranks + 0.5) / n_points
    return design


def cross_2d(origin: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    return float((a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0]))


def convex_hull_indices(points: np.ndarray) -> set[int]:
    if points.shape[0] <= 3:
        return set(range(points.shape[0]))

    order = sorted(range(points.shape[0]), key=lambda i: (points[i, 0], points[i, 1]))
    lower: list[int] = []
    for index in order:
        while len(lower) >= 2 and cross_2d(points[lower[-2]], points[lower[-1]], points[index]) <= 0.0:
            lower.pop()
        lower.append(index)

    upper: list[int] = []
    for index in reversed(order):
        while len(upper) >= 2 and cross_2d(points[upper[-2]], points[upper[-1]], points[index]) <= 0.0:
            upper.pop()
        upper.append(index)

    return set(lower[:-1] + upper[:-1])


def training_indices_without_extrapolation(points: np.ndarray) -> set[int]:
    n_cases = points.shape[0]
    target_train = n_cases - max(1, n_cases // 4)
    train_ids = convex_hull_indices(points)
    for index in range(n_cases):
        if len(train_ids) >= target_train:
            break
        train_ids.add(index)
    return train_ids


def make_parameters(n_cases: int) -> list[dict[str, float | int | bool]]:
    design = latinized_sobol(n_cases)
    alpha = ALPHA_RANGE[0] + design[:, 0] * (ALPHA_RANGE[1] - ALPHA_RANGE[0])
    speed = SPEED_RANGE[0] + design[:, 1] * (SPEED_RANGE[1] - SPEED_RANGE[0])
    points = np.column_stack((alpha, speed))
    train_ids = training_indices_without_extrapolation(points)

    params: list[dict[str, float | int | bool]] = []
    for case_id, (alpha_value, speed_value) in enumerate(points):
        params.append(
            {
                "case": case_id,
                "alpha_deg": float(alpha_value),
                "speed": float(speed_value),
                "train": case_id in train_ids,
            }
        )
    return params


def vector_text(values: np.ndarray) -> str:
    return f"({values[0]:.12g} {values[1]:.12g} {values[2]:.12g})"


def patch_first(pattern: str, repl: str, path: Path) -> None:
    text = path.read_text()
    new_text, count = re.subn(pattern, repl, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"Could not patch {path} with pattern {pattern!r}")
    path.write_text(new_text)


def configure_case(alpha_deg: float, speed: float) -> None:
    alpha = math.radians(alpha_deg)
    drag = np.array([math.cos(alpha), math.sin(alpha), 0.0])
    lift = np.array([-math.sin(alpha), math.cos(alpha), 0.0])
    velocity = speed * drag

    patch_first(
        r"internalField\s+uniform\s+\([^;]+\);",
        f"internalField   uniform {vector_text(velocity)};",
        CASE_DIR / "0.orig" / "U",
    )
    patch_first(
        r"magUInf\s+[-+0-9.eE]+;",
        f"magUInf         {speed:.12g};",
        CASE_DIR / "system" / "controlDict",
    )
    patch_first(
        r"liftDir\s+\([^;]+\);",
        f"liftDir         {vector_text(lift)};",
        CASE_DIR / "system" / "controlDict",
    )
    patch_first(
        r"dragDir\s+\([^;]+\);",
        f"dragDir         {vector_text(drag)};",
        CASE_DIR / "system" / "controlDict",
    )


def processor_count(case_dir: Path) -> int:
    return len([p for p in case_dir.iterdir() if p.is_dir() and p.name.startswith("processor")])


def remove_post_processing(case_dir: Path) -> None:
    shutil.rmtree(case_dir / "postProcessing", ignore_errors=True)
    for processor_dir in case_dir.glob("processor*"):
        shutil.rmtree(processor_dir / "postProcessing", ignore_errors=True)


def numeric_rows(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            rows.append([float(item) for item in stripped.split()])
        except ValueError:
            pass
    if not rows:
        raise RuntimeError(f"No numeric rows in {path}")
    return np.asarray(rows, dtype=float)


def post_files(case_dir: Path, pattern: str) -> list[Path]:
    files = list((case_dir / "postProcessing").glob(pattern))
    for processor_dir in case_dir.glob("processor*"):
        files.extend((processor_dir / "postProcessing").glob(pattern))
    return sorted(files)


def load_surface_files(files: list[Path]) -> np.ndarray:
    if not files:
        return np.empty((0, 0))
    return np.vstack([numeric_rows(path) for path in files])


def save_surface_table(path: Path, surface: np.ndarray) -> None:
    if surface.shape[1] == 4:
        header = "x,y,z,p"
    else:
        header = ",".join(f"col_{i}" for i in range(surface.shape[1]))
    np.savetxt(path, surface, delimiter=",", header=header, comments="")


def time_name_value(path: Path) -> float:
    try:
        return float(path.parent.name)
    except ValueError:
        return -np.inf


def surface_pressure_files_by_time(case_dir: Path) -> list[tuple[float, list[Path]]]:
    files = post_files(case_dir, "surfacePressure/**/*.raw")
    if not files:
        files = post_files(case_dir, "surfaces/**/*.raw")

    groups: dict[float, list[Path]] = {}
    for path in files:
        groups.setdefault(time_name_value(path), []).append(path)

    return [(time, paths) for time, paths in sorted(groups.items())]


def surface_pressure_by_time(case_dir: Path) -> list[tuple[float, np.ndarray]]:
    return [
        (time, load_surface_files(paths))
        for time, paths in surface_pressure_files_by_time(case_dir)
    ]


def latest_force_coefficients(case_dir: Path) -> dict[str, float]:
    files = post_files(case_dir, "forceCoeffs/**/*.dat")
    if not files:
        raise RuntimeError("No force coefficient file was written")
    data = numeric_rows(files[-1])[-1]
    return {
        "Cd": float(data[1]) if data.size > 1 else np.nan,
        "Cl": float(data[3]) if data.size > 3 else np.nan,
        "CmPitch": float(data[5]) if data.size > 5 else np.nan,
    }


def latest_surface_pressure(case_dir: Path) -> np.ndarray:
    surfaces = surface_pressure_by_time(case_dir)
    if not surfaces:
        raise RuntimeError("No raw surface-pressure file was written")
    return surfaces[-1][1]


def copy_latest_surface_pressure_files(case_dir: Path, output_dir: Path) -> None:
    surface_files = surface_pressure_files_by_time(case_dir)
    if not surface_files:
        raise RuntimeError("No raw surface-pressure file was written")
    _, files = surface_files[-1]
    copy_surface_pressure_files(files, case_dir, output_dir)


def copy_surface_pressure_files(files: list[Path], case_dir: Path, output_dir: Path) -> None:
    for path in files:
        destination = output_dir / path.relative_to(case_dir)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)


def archive_run(case_id: int, parameters: dict[str, float | int | bool]) -> dict[str, float]:
    run_dir = RESULTS_DIR / "cases" / f"run_{case_id:03d}"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    shutil.copytree(CASE_DIR / "exported_data", run_dir / "exported_data")

    coeffs = latest_force_coefficients(CASE_DIR)
    surface = latest_surface_pressure(CASE_DIR)
    np.save(run_dir / "surface_pressure.npy", surface)
    save_surface_table(run_dir / "surface_pressure.csv", surface)
    copy_latest_surface_pressure_files(CASE_DIR, run_dir / "surface_pressure_raw")

    with (run_dir / "parameters.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case", "alpha_deg", "speed", "train"])
        writer.writeheader()
        writer.writerow(parameters)

    return coeffs


def processor_files(data_dir: Path, field: str) -> list[tuple[int, Path]]:
    pattern = re.compile(rf"^{re.escape(field)}_proc_(\d+)\.npy$")
    files = []
    for path in (data_dir / field).glob(f"{field}_proc_*.npy"):
        match = pattern.match(path.name)
        if match:
            files.append((int(match.group(1)), path))
    files.sort(key=lambda item: item[0])
    if not files:
        raise RuntimeError(f"No {field} processor files in {data_dir}")
    return files


def load_snapshot(run_dir: Path) -> tuple[list[int], list[int], dict[str, np.ndarray], np.ndarray]:
    data_dir = run_dir / "exported_data"
    p_parts = [(proc, np.load(path)) for proc, path in processor_files(data_dir, "p")]
    u_parts = [(proc, np.load(path)) for proc, path in processor_files(data_dir, "U")]
    v_parts = [(proc, np.load(path)) for proc, path in processor_files(data_dir, "cellVolumes")]
    proc_ids = [proc for proc, _ in p_parts]
    if (
        proc_ids != [proc for proc, _ in u_parts]
        or proc_ids != [proc for proc, _ in v_parts]
    ):
        raise RuntimeError("Processor ids differ between exported fields")

    counts = [part.shape[0] for _, part in p_parts]
    p = np.concatenate([part[:, 0] for _, part in p_parts])
    u = np.concatenate([part[:, :, 0] for _, part in u_parts])
    volumes = np.concatenate([part[:, 0] for _, part in v_parts])
    fields = {"p": p, "Ux": u[:, 0], "Uy": u[:, 1]}
    return proc_ids, counts, fields, volumes


def stack_raw_state(fields: dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate((fields["p"], fields["Ux"], fields["Uy"]))


def state_reference_vector(n_cells: int) -> np.ndarray:
    return np.concatenate(
        (
            np.full(n_cells, REF_PRESSURE),
            np.full(n_cells, REF_SPEED),
            np.full(n_cells, REF_SPEED),
        )
    )


def minmax_scale_state(
    raw_state: np.ndarray,
    train_ids: list[int],
    n_cells: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    reference = state_reference_vector(n_cells)[:, np.newaxis]
    nondimensional_state = raw_state / reference
    train_state = nondimensional_state[:, train_ids]

    data_min = np.zeros((3 * n_cells, 1), dtype=float)
    data_range = np.ones((3 * n_cells, 1), dtype=float)
    for start, stop in (
        (0, n_cells),
        (n_cells, 2 * n_cells),
        (2 * n_cells, 3 * n_cells),
    ):
        block = train_state[start:stop]
        block_min = float(np.min(block))
        block_range = float(np.max(block) - block_min)
        if block_range == 0.0:
            block_range = 1.0
        data_min[start:stop] = block_min
        data_range[start:stop] = block_range

    data_range[data_range == 0.0] = 1.0
    scaled_state = (nondimensional_state - data_min) / data_range
    return scaled_state, data_min, data_range, reference


def inverse_scale_state(
    scaled_state: np.ndarray,
    data_min: np.ndarray,
    data_range: np.ndarray,
    reference: np.ndarray,
) -> np.ndarray:
    return (scaled_state * data_range + data_min) * reference


def normalize_parameters(parameters: list[dict[str, float | int | bool]]) -> np.ndarray:
    x = np.empty((len(parameters), 2), dtype=float)
    for row, item in enumerate(parameters):
        x[row, 0] = 2 * (float(item["alpha_deg"]) - ALPHA_RANGE[0]) / (ALPHA_RANGE[1] - ALPHA_RANGE[0]) - 1
        x[row, 1] = 2 * (float(item["speed"]) - SPEED_RANGE[0]) / (SPEED_RANGE[1] - SPEED_RANGE[0]) - 1
    return x


def thin_plate_kernel(r: np.ndarray) -> np.ndarray:
    kernel = np.zeros_like(r)
    positive = r > 0.0
    kernel[positive] = (r[positive] ** 2) * np.log(r[positive])
    return kernel


def pairwise_distances(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.linalg.norm(x[:, np.newaxis, :] - y[np.newaxis, :, :], axis=2)


def fit_thin_plate_spline(
    x_train: np.ndarray,
    y_train: np.ndarray,
    smoothing: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_samples = x_train.shape[0]
    kernel = thin_plate_kernel(pairwise_distances(x_train, x_train))
    kernel += max(0.0, smoothing) * np.eye(n_samples)
    affine_basis = np.column_stack((np.ones(n_samples), x_train))
    system = np.block(
        [
            [kernel, affine_basis],
            [affine_basis.T, np.zeros((affine_basis.shape[1], affine_basis.shape[1]))],
        ]
    )
    rhs = np.vstack((y_train, np.zeros((affine_basis.shape[1], y_train.shape[1]))))
    solution, *_ = np.linalg.lstsq(system, rhs, rcond=None)
    weights = solution[:n_samples]
    affine = solution[n_samples:]
    return x_train.copy(), weights, affine


def evaluate_thin_plate_spline(model: tuple[np.ndarray, np.ndarray, np.ndarray], x_eval: np.ndarray) -> np.ndarray:
    x_train, weights, affine = model
    kernel = thin_plate_kernel(pairwise_distances(x_eval, x_train))
    affine_basis = np.column_stack((np.ones(x_eval.shape[0]), x_eval))
    return kernel @ weights + affine_basis @ affine


def write_field_parts(
    predictions: np.ndarray,
    proc_ids: list[int],
    counts: list[int],
    output_dir: Path,
) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    n_cells = sum(counts)
    p = predictions[:n_cells]
    ux = predictions[n_cells : 2 * n_cells]
    uy = predictions[2 * n_cells :]
    (output_dir / "p").mkdir(parents=True)
    (output_dir / "U").mkdir(parents=True)

    start = 0
    for proc, count in zip(proc_ids, counts):
        stop = start + count
        np.save(output_dir / "p" / f"p_proc_{proc}.npy", np.asfortranarray(p[start:stop]))
        u_part = np.zeros((count, 3, predictions.shape[1]), dtype=float, order="F")
        u_part[:, 0, :] = ux[start:stop]
        u_part[:, 1, :] = uy[start:stop]
        np.save(output_dir / "U" / f"U_proc_{proc}.npy", u_part)
        start = stop


def update_prediction_dict(n_times: int) -> None:
    text = (PODI_CASE / "system" / "numpyToFoamDict").read_text()
    text = re.sub(r"endTime\s+[-+0-9.eE]+;", f"endTime     {n_times};", text, count=1)
    (PODI_CASE / "system" / "numpyToFoamDict").write_text(text)


def time_selector(n_times: int) -> str:
    return ",".join(str(i) for i in range(1, n_times + 1))


def copy_auxiliary_field_to_prediction_times(case_dir: Path, field: str, n_times: int) -> None:
    for processor_dir in case_dir.glob("processor*"):
        source = processor_dir / "0" / field
        if not source.exists():
            raise FileNotFoundError(f"Required auxiliary field is missing: {source}")

        source_text = source.read_text()
        for time_index in range(1, n_times + 1):
            destination = processor_dir / str(time_index) / field
            destination.parent.mkdir(parents=True, exist_ok=True)
            text = re.sub(
                r'location\s+"[^"]+";',
                f'location    "{time_index}";',
                source_text,
                count=1,
            )
            destination.write_text(text)


def ignore_podi_case_files(directory: str, names: list[str]) -> set[str]:
    # Keep system/ dictionaries intact; skip generated fields and logs only.
    return {
        name
        for name in names
        if name in {"500", "postProcessing", "exported_data"}
        or name == "phi"
        or name.startswith("log.")
    }


def prepare_podi_case(source: Path, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination, ignore=ignore_podi_case_files)
    remove_post_processing(destination)


def write_podi_surface_dict(case_dir: Path) -> None:
    text = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      podiSurfacePressureDict;
}

functions
{
    surfacePressure
    {
        type            surfaces;
        libs            (sampling);
        surfaceFormat   raw;
        interpolationScheme cellPoint;
        fields          (p);
        executeControl  timeStep;
        writeControl    timeStep;
        surfaces
        (
            airfoil
            {
                type    patch;
                patches (walls);
            }
        );
    }
}
"""
    (case_dir / "system" / "podiSurfacePressureDict").write_text(text)


def predicted_surface_pressures(
    evaluation_ids: list[int],
    surfaces_by_time: list[tuple[float, np.ndarray]],
) -> dict[int, np.ndarray]:
    surfaces = {}
    for local_index, case_id in enumerate(evaluation_ids):
        if local_index < len(surfaces_by_time):
            surfaces[case_id] = surfaces_by_time[local_index][1]
        else:
            surfaces[case_id] = np.empty((0, 0))
    return surfaces


def align_surface_pressure(truth: np.ndarray, prediction: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if truth.shape[1] < 4 or prediction.shape[1] < 4:
        n = min(truth.shape[0], prediction.shape[0])
        coords = np.column_stack((np.arange(n, dtype=float), np.zeros(n)))
        return coords, truth[:n, -1], prediction[:n, -1]

    truth_keys = np.round(truth[:, :3], decimals=10)
    prediction_keys = np.round(prediction[:, :3], decimals=10)
    prediction_lookup = {tuple(key): i for i, key in enumerate(prediction_keys)}

    truth_indices = []
    prediction_indices = []
    for i, key in enumerate(truth_keys):
        j = prediction_lookup.get(tuple(key))
        if j is not None:
            truth_indices.append(i)
            prediction_indices.append(j)

    if not truth_indices:
        n = min(truth.shape[0], prediction.shape[0])
        truth_indices = list(range(n))
        prediction_indices = list(range(n))

    truth_idx = np.asarray(truth_indices, dtype=int)
    prediction_idx = np.asarray(prediction_indices, dtype=int)
    return truth[truth_idx, :2], truth[truth_idx, -1], prediction[prediction_idx, -1]


def surface_pressure_coefficient_table(
    case_id: int,
    parameters: list[dict[str, float | int | bool]],
    prediction: np.ndarray,
) -> np.ndarray:
    truth = np.load(RESULTS_DIR / "cases" / f"run_{case_id:03d}" / "surface_pressure.npy")
    coords, truth_p, prediction_p = align_surface_pressure(truth, prediction)
    x = coords[:, 0]
    x_min = float(np.min(x))
    x_range = float(np.max(x) - x_min)
    if x_range == 0.0:
        x_range = 1.0
    x_normalized = (x - x_min) / x_range
    dynamic_pressure = 0.5 * float(parameters[case_id]["speed"]) ** 2
    truth_negative_cp = -(truth_p / dynamic_pressure)
    prediction_negative_cp = -(prediction_p / dynamic_pressure)
    difference = prediction_negative_cp - truth_negative_cp
    return np.column_stack((x_normalized, truth_negative_cp, prediction_negative_cp, difference))


def save_surface_pressure_coefficient_table(path: Path, table: np.ndarray) -> None:
    header = "x_normalized,minus_Cp_openfoam,minus_Cp_podi,difference"
    np.savetxt(path, table, delimiter=",", header=header, comments="")


def plot_surface_pressure_comparisons(
    evaluation_ids: list[int],
    parameters: list[dict[str, float | int | bool]],
    predicted_surfaces: dict[int, np.ndarray],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir = RESULTS_DIR / "surface_pressure_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for case_id in evaluation_ids:
        prediction = predicted_surfaces[case_id]
        if prediction.size == 0:
            print(f"No predicted surface pressure found for case {case_id}")
            continue

        cp_table = surface_pressure_coefficient_table(case_id, parameters, prediction)
        x_normalized = cp_table[:, 0]
        order = np.argsort(x_normalized)
        truth_negative_cp = cp_table[:, 1]
        prediction_negative_cp = cp_table[:, 2]
        difference = cp_table[:, 3]

        fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.2), sharex=True, constrained_layout=True)
        split = "training" if parameters[case_id]["train"] else "test"
        title = (
            f"case {case_id:03d} ({split}): "
            f"alpha={float(parameters[case_id]['alpha_deg']):.2f} deg, "
            f"U={float(parameters[case_id]['speed']):.2f}"
        )
        fig.suptitle(title)

        axes[0].plot(x_normalized[order], truth_negative_cp[order], "k.", markersize=4, label="OpenFOAM")
        axes[0].plot(x_normalized[order], prediction_negative_cp[order], "C3.", markersize=4, label="PODI")
        axes[0].set_ylabel("-Cp")
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()

        axes[1].plot(x_normalized[order], difference[order], "C0.", markersize=4)
        axes[1].axhline(0.0, color="k", linewidth=1.0, alpha=0.6)
        axes[1].set_xlabel("normalized x")
        axes[1].set_ylabel("PODI - OpenFOAM")
        axes[1].grid(True, alpha=0.3)
        fig.savefig(plot_dir / f"surface_pressure_case_{case_id:03d}.png", dpi=160)
        plt.close(fig)


def plot_force_coefficient_comparisons(rows: list[dict[str, float | int | bool]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir = RESULTS_DIR / "force_coefficient_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    cases = np.array([row["case"] for row in rows], dtype=int)
    is_train = np.array([row["train"] for row in rows], dtype=bool)
    fields = (
        ("Cd", "drag coefficient"),
        ("Cl", "lift coefficient"),
        ("CmPitch", "pitching moment coefficient"),
    )

    fig, axes = plt.subplots(len(fields), 1, figsize=(7.0, 6.2), sharex=True, constrained_layout=True)
    for ax, (name, label) in zip(axes, fields):
        truth = np.array([row[f"{name}_true"] for row in rows], dtype=float)
        podi = np.array([row[f"{name}_podi"] for row in rows], dtype=float)
        ax.plot(cases, truth, "k-", linewidth=1.0, alpha=0.4)
        ax.plot(cases, podi, "C3--", linewidth=1.0, alpha=0.4)
        ax.scatter(cases[is_train], truth[is_train], c="k", marker="o", label="OpenFOAM train")
        ax.scatter(cases[~is_train], truth[~is_train], c="k", marker="x", label="OpenFOAM test")
        ax.scatter(cases[is_train], podi[is_train], c="C3", marker="s", label="PODI train")
        ax.scatter(cases[~is_train], podi[~is_train], c="C3", marker="+", label="PODI test")
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncols=2)
    axes[-1].set_xlabel("case")
    fig.savefig(plot_dir / "force_coefficients_by_case.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(fields), figsize=(9.0, 3.2), constrained_layout=True)
    for ax, (name, label) in zip(axes, fields):
        truth = np.array([row[f"{name}_true"] for row in rows], dtype=float)
        podi = np.array([row[f"{name}_podi"] for row in rows], dtype=float)
        lower = float(min(np.min(truth), np.min(podi)))
        upper = float(max(np.max(truth), np.max(podi)))
        if lower == upper:
            lower -= 1.0
            upper += 1.0
        ax.plot([lower, upper], [lower, upper], "k:", linewidth=1.0)
        ax.scatter(truth[is_train], podi[is_train], c=cases[is_train], cmap="viridis", s=36, marker="s", label="train")
        ax.scatter(truth[~is_train], podi[~is_train], c=cases[~is_train], cmap="viridis", s=44, marker="x", label="test")
        ax.set_xlabel(f"OpenFOAM {name}")
        ax.set_ylabel(f"PODI {name}")
        ax.set_title(label, fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.savefig(plot_dir / "force_coefficients_parity.png", dpi=160)
    plt.close(fig)


def train_and_compare(
    parameters: list[dict[str, float | int | bool]],
    coeff_rows: list[dict[str, float]],
    spline_smoothing: float,
    modes: int,
) -> None:
    run_dirs = [RESULTS_DIR / "cases" / f"run_{i:03d}" for i in range(len(parameters))]
    snapshots = []
    proc_ids = []
    counts = []
    volumes = None
    for run_dir in run_dirs:
        proc_ids, counts, fields, volumes = load_snapshot(run_dir)
        snapshots.append(stack_raw_state(fields))
    raw_state_matrix = np.column_stack(snapshots)
    assert volumes is not None

    train_ids = [int(item["case"]) for item in parameters if item["train"]]
    evaluation_ids = list(range(len(parameters)))
    x_all = normalize_parameters(parameters)
    x_train = x_all[train_ids]
    x_evaluation = x_all[evaluation_ids]

    n_cells = sum(counts)
    scaled_state_matrix, data_min, data_range, reference = minmax_scale_state(
        raw_state_matrix,
        train_ids,
        n_cells,
    )
    mean_state = scaled_state_matrix[:, train_ids].mean(axis=1, keepdims=True)
    centered = scaled_state_matrix[:, train_ids] - mean_state
    weights = np.sqrt(np.tile(volumes, 3))[:, np.newaxis]
    weighted_centered = weights * centered
    weighted_modes, singular_values, vt = np.linalg.svd(weighted_centered, full_matrices=False)
    max_modes = weighted_modes.shape[1]
    n_modes = max_modes if modes <= 0 else min(modes, max_modes)
    pod_modes = weighted_modes[:, :n_modes] / weights
    coefficients = np.diag(singular_values[:n_modes]) @ vt[:n_modes, :]

    coefficient_spline = fit_thin_plate_spline(x_train, coefficients.T, spline_smoothing)
    predicted_coeffs = evaluate_thin_plate_spline(coefficient_spline, x_evaluation)
    predicted_scaled_states = mean_state + pod_modes @ predicted_coeffs.T
    predicted_states = inverse_scale_state(predicted_scaled_states, data_min, data_range, reference)

    force_targets = np.array([[row["Cd"], row["Cl"], row["CmPitch"]] for row in coeff_rows], dtype=float)
    force_spline = fit_thin_plate_spline(x_train, force_targets[train_ids], spline_smoothing)
    force_predictions = evaluate_thin_plate_spline(force_spline, x_evaluation)

    write_field_parts(predicted_states, proc_ids, counts, PREDICTION_DATA)
    prepare_podi_case(CASE_DIR, PODI_CASE)
    shutil.copytree(PREDICTION_DATA, PODI_CASE / "podi_predictions")
    update_prediction_dict(len(evaluation_ids))
    write_podi_surface_dict(PODI_CASE)
    nproc = processor_count(PODI_CASE)
    run(["mpirun", "-np", str(nproc), "numpyToFoam", "-parallel"], cwd=PODI_CASE)
    copy_auxiliary_field_to_prediction_times(PODI_CASE, "nut", len(evaluation_ids))
    copy_auxiliary_field_to_prediction_times(PODI_CASE, "nuTilda", len(evaluation_ids))
    run(
        [
            "mpirun",
            "-np",
            str(nproc),
            "simpleFoam",
            "-postProcess",
            "-parallel",
            "-dict",
            "system/podiSurfacePressureDict",
            "-time",
            time_selector(len(evaluation_ids)),
        ],
        cwd=PODI_CASE,
    )

    predicted_surface_files = surface_pressure_files_by_time(PODI_CASE)
    if not predicted_surface_files:
        raise RuntimeError(
            "PODI surface-pressure sampling wrote no raw files; "
            "check the simpleFoam -postProcess output for surfacePressure"
    )
    predicted_surfaces = predicted_surface_pressures(
        evaluation_ids,
        [(time, load_surface_files(paths)) for time, paths in predicted_surface_files],
    )
    for case_id, surface in predicted_surfaces.items():
        run_dir = RESULTS_DIR / "cases" / f"run_{case_id:03d}"
        if surface.size:
            np.save(run_dir / "surface_pressure_podi.npy", surface)
            save_surface_table(run_dir / "surface_pressure_podi.csv", surface)
            cp_table = surface_pressure_coefficient_table(case_id, parameters, surface)
            save_surface_pressure_coefficient_table(
                run_dir / "surface_pressure_coefficient_comparison.csv",
                cp_table,
            )
    for local_index, case_id in enumerate(evaluation_ids):
        if local_index < len(predicted_surface_files):
            _, files = predicted_surface_files[local_index]
            copy_surface_pressure_files(
                files,
                PODI_CASE,
                RESULTS_DIR / "cases" / f"run_{case_id:03d}" / "surface_pressure_podi_raw",
            )
    plot_surface_pressure_comparisons(evaluation_ids, parameters, predicted_surfaces)
    rows = []
    for local_index, case_id in enumerate(evaluation_ids):
        true_state = raw_state_matrix[:, case_id]
        pred_state = predicted_states[:, local_index]
        rows.append(
            {
                "case": case_id,
                "train": parameters[case_id]["train"],
                "alpha_deg": parameters[case_id]["alpha_deg"],
                "speed": parameters[case_id]["speed"],
                "Cd_true": coeff_rows[case_id]["Cd"],
                "Cd_podi": force_predictions[local_index, 0],
                "Cl_true": coeff_rows[case_id]["Cl"],
                "Cl_podi": force_predictions[local_index, 1],
                "CmPitch_true": coeff_rows[case_id]["CmPitch"],
                "CmPitch_podi": force_predictions[local_index, 2],
                "field_relative_l2": np.linalg.norm(true_state - pred_state) / np.linalg.norm(true_state),
            }
        )

    RESULTS_DIR.mkdir(exist_ok=True)
    with (RESULTS_DIR / "comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    plot_force_coefficient_comparisons(rows)
    np.save(RESULTS_DIR / "singular_values.npy", singular_values)
    with (RESULTS_DIR / "scaling_reference.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reference_speed", "reference_pressure"])
        writer.writeheader()
        writer.writerow(
            {
                "reference_speed": REF_SPEED,
                "reference_pressure": REF_PRESSURE,
            }
        )
    print(f"POD modes used: {n_modes}")
    print(f"Thin-plate spline smoothing: {spline_smoothing:g}")
    print(f"Wrote {RESULTS_DIR / 'comparison.csv'}")
    print(f"Wrote {RESULTS_DIR / 'surface_pressure_plots'}")
    print(f"Wrote {RESULTS_DIR / 'force_coefficient_plots'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=8, help="Number of parameter cases")
    parser.add_argument("--spline-smoothing", type=float, default=0.0, help="Thin-plate spline smoothing")
    parser.add_argument("--degree", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--modes", type=int, default=0, help="POD modes; 0 keeps all train modes")
    args = parser.parse_args()

    if args.cases < MIN_CASES:
        raise SystemExit(f"--cases must be at least {MIN_CASES}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    parameters = make_parameters(args.cases)
    coeff_rows = []

    with (RESULTS_DIR / "parameters.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case", "alpha_deg", "speed", "train"])
        writer.writeheader()
        writer.writerows(parameters)

    for item in parameters:
        case_id = int(item["case"])
        print(f"\n=== Case {case_id:03d}: alpha={item['alpha_deg']:.3f} deg, speed={item['speed']:.3f} ===")
        configure_case(float(item["alpha_deg"]), float(item["speed"]))
        run(["./Allclean"], cwd=CASE_DIR)
        run(["./Allrun"], cwd=CASE_DIR)
        nproc = processor_count(CASE_DIR)
        if nproc == 0:
            raise RuntimeError("The airfoil case did not create processor directories")
        shutil.rmtree(CASE_DIR / "exported_data", ignore_errors=True)
        run(["mpirun", "-np", str(nproc), "foamToNumpy", "-parallel"], cwd=CASE_DIR)
        coeff_rows.append(archive_run(case_id, item))

    train_and_compare(parameters, coeff_rows, args.spline_smoothing, args.modes)


if __name__ == "__main__":
    main()
