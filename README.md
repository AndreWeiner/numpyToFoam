# numpyToFoam

## General Description
`numpyToFoam` is a custom OpenFOAM utility for reconstructing OpenFOAM field data from NumPy (`.npy`) arrays. It supports scalar, vector, symmetric tensor, and full tensor fields.

The field values over the computational cells at a given time instant constitute a **snapshot**. A collection of such snapshots over multiple time steps forms the input dataset used by the utility.

These snapshots are processed sequentially in time, avoiding the need to store large datas

---

## Dependencies
- OpenFOAM 2112, 2206, 2312, 2412, or 2512
- A working C++ compiler available in the OpenFOAM environment (e.g. `g++`)

---

## Compilation and Installation
In a working OpenFOAM environment:

```bash
./Allwmake
```

Both utilities will be compiled and made available within your OpenFOAM environment.

To clean both builds:

```bash
./Allwclean
```

To compile only `numpyToFoam`:

```bash
cd src/numpyToFoam
wmake
```

To clean only `numpyToFoam`:

```bash
cd src/numpyToFoam
wclean
```

---

## Usage

### Requirements

To use the utility, the following are required:

#### 1. OpenFOAM Case
- A decomposed OpenFOAM case (parallel setup required)
- Mesh must be available in `processor*/` directories
- Optional: `0/` directory with initial and boundary conditions

#### 2. NumPy Snapshot Files
- Data must be stored as `.npy` files, organised in one subdirectory per field:

```
data/
├── p/
│   ├── p_proc_0.npy
│   └── p_proc_1.npy
└── U/
    ├── U_proc_0.npy
    └── U_proc_1.npy
```

- Naming convention:

```bash
dataDir/fieldName/fieldName_proc_i.npy
```

where:
- `fieldName` is the OpenFOAM field name, for example `p` or `U`
- `i` is the processor index

#### Supported Shapes

| Field Type        | Shape                     |
|-------------------|---------------------------|
| Scalar            | `(nCells, nTimeSteps)`    |
| Vector            | `(nCells, 3, nTimeSteps)` |
| Symmetric Tensor  | `(nCells, 6, nTimeSteps)` |
| Tensor            | `(nCells, 9, nTimeSteps)` |

#### Storage Format
- Both **row-major (C-order)** and **column-major (Fortran-order)** storage formats are supported
- The utility automatically detects the array storage format from the `.npy` header
- **Recommended:** use column-major format for faster snapshot access

#### 3. `numpyToFoamDict`
A `numpyToFoamDict` file must be present in the `system/` directory.

This dictionary controls:
- the location of the `.npy` files (`dataDir`)
- the list of fields to reconstruct
- the time range and time-step size used for writing OpenFOAM fields

Example:

```plaintext
dataDir    data;

fields     (p U);

time
{
    startTime   0;
    endTime     1;
    deltaT      0.1;
}
```

#### Dictionary entries

| Entry | Description |
|-------|-------------|
| `dataDir` | Directory containing the `.npy` files, relative to the case root (or absolute). Default: `data`. |
| `fields` | List of OpenFOAM field names to reconstruct (e.g. `(p U)`). The storage format and precision are detected automatically from each `.npy` header. |
| `time/startTime` | Time value of the first snapshot to write. |
| `time/endTime` | Time value of the last snapshot to write. |
| `time/deltaT` | Time-step size used to generate uniformly spaced output times between `startTime` and `endTime`. |

The number of generated time steps must match the number of snapshots stored in the `.npy` files.

---

### Running the Utility

```bash
mpirun -np 2 numpyToFoam -parallel
```

Replace `2` with the number of processors used for case decomposition.

---

# foamToNumpy

## General Description
`foamToNumpy` is a custom OpenFOAM utility for extracting OpenFOAM field data into NumPy (`.npy`) arrays. It is the reverse of `numpyToFoam`. It supports scalar, vector, symmetric tensor, and full tensor fields.

All snapshots for a given field are packed into a single `.npy` file per processor, avoiding repeated file I/O.

The utility can optionally export mesh geometry data (cell centres, cell volumes) and the selected time values as additional `.npy` files.

---

## Dependencies
- OpenFOAM 2112, 2206, 2312, 2412, or 2512
- A working C++ compiler available in the OpenFOAM environment (e.g. `g++`)

---

## Compilation and Installation
In a working OpenFOAM environment:

```bash
./Allwmake
```

Both utilities will be compiled and made available within your OpenFOAM environment.

To clean both builds:

```bash
./Allwclean
```

To compile only `foamToNumpy`:

```bash
cd src/foamToNumpy
wmake
```

To clean only `foamToNumpy`:

```bash
cd src/foamToNumpy
wclean
```

---

## Usage

### Requirements

#### 1. OpenFOAM Case
- A decomposed OpenFOAM case (parallel setup required)
- Field data must be available in `processor*/` directories for the selected time range

#### 2. Output File Structure
Each field is written into its own subdirectory inside `dataDir`:

```
data/
├── p/
│   ├── p_proc_0.npy
│   └── p_proc_1.npy
└── U/
    ├── U_proc_0.npy
    └── U_proc_1.npy
```

Naming convention:

```bash
dataDir/fieldName/fieldName_proc_i.npy
```

where `i` is the processor index.

#### Output Shapes

| Field Type        | Shape                      |
|-------------------|----------------------------|
| Scalar            | `(nCells, nTimeSteps)`     |
| Vector            | `(nCells, 3, nTimeSteps)`  |
| Symmetric Tensor  | `(nCells, 6, nTimeSteps)`  |
| Tensor            | `(nCells, 9, nTimeSteps)`  |

#### Storage Format
- Both **row-major (C-order)** and **column-major (Fortran-order)** storage formats are supported
- The storage format is set via `storageOrder` in the dictionary
- **Recommended:** use column-major format for faster snapshot access along the time axis

#### 3. `foamToNumpyDict`
A `foamToNumpyDict` file must be present in the `system/` directory.

This dictionary controls:
- the output location for `.npy` files (`dataDir`)
- the list of fields to extract and their output data type
- the time range and sub-sampling stride
- optional geometry and time exports (`exportData`)
- the array storage order

Example:

```plaintext
dataDir       data;

fields
{
    names       (p U);
    dataType    float64;
}

exportData
{
    cellCentre   true;
    cellVolumes  true;
    writeTimes   true;
    dataType     float64;
}

storageOrder   F;         // F or C

time
{
    startTime   0.0;
    endTime     0.5;
    every       1;
}
```

#### Dictionary entries

| Entry | Description |
|-------|-------------|
| `dataDir` | Output directory for `.npy` files, relative to the case root (or absolute). Created automatically if it does not exist. |
| `fields/names` | List of OpenFOAM field names to extract. |
| `fields/dataType` | Numeric type for field output: `float32` or `float64` (default: `float64`). |
| `exportData/cellCentre` | If `true`, exports cell-centre coordinates to `dataDir/cellCentre/cellCentre_proc_i.npy` with shape `(nCells, 3, 1)`. |
| `exportData/cellVolumes` | If `true`, exports cell volumes to `dataDir/cellVolumes/cellVolumes_proc_i.npy` with shape `(nCells, 1)`. |
| `exportData/writeTimes` | If `true`, exports the selected time values to `dataDir/times/times.npy` with shape `(nTimeSteps, 1)`. |
| `exportData/dataType` | Numeric type for geometry/time exports: `float32` or `float64` (default: `float64`). |
| `storageOrder` | Array storage order: `F` (Fortran/column-major) or `C` (C/row-major). |
| `time/startTime` | First time step to include. |
| `time/endTime` | Last time step to include. |
| `time/every` | Sub-sampling stride: write every N-th time step within the range (e.g. `1` keeps all, `2` keeps every other). |

---

> **Note — using both utilities in tandem.**
> When `foamToNumpy` and `numpyToFoam` are used together, the time settings in their respective dictionaries must be consistent: the number of snapshots written by `foamToNumpy` must equal the number of time steps generated by `numpyToFoam`.
>
> Given a case with time steps `0.1, 0.2, 0.3, 0.4, 0.5`, the following pair of settings selects all five snapshots and reconstructs them at the same time labels:
>
> ```plaintext
> # foamToNumpyDict          # numpyToFoamDict
> time                       time
> {                          {
>     startTime   0.1;           startTime   0.1;
>     endTime     0.5;           endTime     0.5;
>     every       1;             deltaT      0.1;
> }                          }
> ```
>
> If `every 2` is used in `foamToNumpyDict` (selecting `0.1, 0.3, 0.5` — three snapshots), then `numpyToFoamDict` must generate exactly three time steps, e.g. `startTime 0.1`, `endTime 0.5`, `deltaT 0.2`.

---

### Running the Utility

```bash
mpirun -np 2 foamToNumpy -parallel
```

Replace `2` with the number of processors used for case decomposition.

> **Note — skip-on-conflict behaviour**
> If the output directory for a field or an `exportData` entry already exists, that export is skipped with a warning and no data is overwritten. Remove or rename the existing output directories before re-running to avoid silent skips.
---

# For Developers

A unit test framework is provided to validate both utilities across multiple OpenFOAM versions using Apptainer containers. The pipeline runs an `icoFoam` cavity simulation, exports the results to `.npy` with `foamToNumpy`, reconstructs the fields with `numpyToFoam`, and verifies correctness by comparing MD5 checksums of the reconstructed processor field files against the original simulation output.

**Tested OpenFOAM versions:** 2112, 2206, 2312, 2412, 2512

For each version the following steps are executed in order:

| Step | Description |
|------|-------------|
| `foamToNumpy build` | Compiles `foamToNumpy` inside the container |
| `numpyToFoam build` | Compiles `numpyToFoam` inside the container |
| `Allrun` | Runs the `icoFoam` cavity simulation |
| `foamToNumpy run` | Exports simulation fields to `.npy` |
| `Clean_proc_data` | Removes the original OpenFOAM field data from the processor directories |
| `numpyToFoam run` | Reconstructs fields from `.npy` |
| `checksum match` | Compares MD5 checksums of reconstructed vs. original fields |

A pass/fail status is printed for each step and version at the end of the run.

---

## Local

### Prerequisites
- [Apptainer](https://apptainer.org) installed

### Running

```bash
cd unittest
bash unittest.sh
```

On the **first run**, the script builds an Apptainer `.sif` image for each OpenFOAM version from Docker Hub and caches them in `unittest/of_versions/`. Subsequent runs reuse the cached images and skip the build step. Logs for each step are written to `unittest/run/of{version}/`.

---

## GitHub Actions

The CI pipeline is defined in [`.github/workflows/main.yml`](.github/workflows/main.yml) and triggers automatically on every push and pull request. It runs `unittest/actions_unittest.sh`, which pulls pre-built Apptainer images from the GitHub Container Registry (GHCR) instead of building them locally, making it faster.

To monitor a run, go to the **Actions** tab of the repository, select the `numpyToFoamTest` workflow, and open the latest run. If any step fails, the last few lines of every relevant log file are printed and the workflow exits with a non-zero status.

---



# Limitations

1. **Parallel-only execution** (`numpyToFoam`, `foamToNumpy`)
   - Serial execution is not supported by either utility

2. **Precision — `numpyToFoam`**
   - If the input `.npy` files are stored in single precision (`float32`), the values are converted to OpenFOAM `scalar` precision during reading
   - In most OpenFOAM builds, this means conversion to **double precision**

3. **No in-built post-processing** (`numpyToFoam`)
   - Derived quantities cannot be computed directly within the utility
   - Post-processing must be performed separately after reconstruction

4. **No missing-time handling** (`numpyToFoam`, `foamToNumpy`)
   - Data must exist for every time step defined in the `time` settings. If any time step has missing data, the utility will terminate abruptly without graceful exception handling

5. **No region-based data write** (`numpyToFoam`, `foamToNumpy`)
   - Writing data for specific mesh regions is not currently supported
   - This feature is planned for a future update

6. **Boundary data not supported** (`numpyToFoam`, `foamToNumpy`)
   - Only internal cell data is read and written; boundary patch field values are not handled

7. **Finite area fields not supported** (`numpyToFoam`, `foamToNumpy`)
   - Only volumetric field types are supported;
