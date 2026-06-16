# Tutorial

## Mode Visualization and Reconstruction Using `foamToNumpy` and `numpyToFoam`

One useful application of these utilities is the visualization of POD modes computed from a snapshot matrix, as well as the visualization of low-order reconstructed fields. In this tutorial, `foamToNumpy` first exports pressure snapshots from an OpenFOAM cavity-flow simulation. A small NumPy script then computes POD modes and a low-rank reconstruction. Finally, `numpyToFoam` writes those arrays back into OpenFOAM format for visualization.

The tutorial is based on the flow-across-cavity case provided in `tutorials/of_cavity`.

You can run the complete workflow with:

```bash
cd tutorials
./run_tutorial.sh
```

To remove generated tutorial data and cases:

```bash
./reset_tutorial.sh
```

The script builds both utilities, runs the cavity simulation, exports pressure and velocity snapshots with `foamToNumpy`, computes the POD data with NumPy, and creates two output cases:

- `tutorials/of_cavity_reconstructed`
- `tutorials/of_cavity_modes`

---

## 1. Run the OpenFOAM Simulation

In a working OpenFOAM environment, run:

```bash
cd tutorials/of_cavity
./Allrun
```

This runs the simulation and stores decomposed field data in the `processor*` directories. The helper script detects the number of processor directories automatically after this step.

---

## 2. Export Snapshot Data with `foamToNumpy`

The tutorial case contains `system/foamToNumpyDict`, which exports the pressure field `p` and velocity field `U` from time `0.01` to `0.5`:

```plaintext
dataDir       exported_data;

fields
{
    names       (p U);
    dataType    float64;
}

exportData
{
    cellCentre   false;
    cellVolumes  true;
    writeTimes   true;
    dataType     float64;
}

storageOrder   F;

time
{
    startTime   0.01;
    endTime     0.5;
    every       1;
}
```

Run `foamToNumpy` in parallel from the case directory:

```bash
cd tutorials/of_cavity
NPROC="$(find . -maxdepth 1 -type d -name 'processor*' | wc -l)"
mpirun -np "$NPROC" foamToNumpy -parallel
```

This creates:

```text
exported_data/
├── p/
│   ├── p_proc_0.npy
│   ├── p_proc_1.npy
│   └── ...
├── U/
│   ├── U_proc_0.npy
│   ├── U_proc_1.npy
│   └── ...
├── cellVolumes/
│   ├── cellVolumes_proc_0.npy
│   ├── cellVolumes_proc_1.npy
│   └── ...
└── times/
    └── times.npy
```

Each `p_proc_i.npy` file contains the pressure snapshots for one processor partition with shape `(nCells, nTimes)`. Each `U_proc_i.npy` file contains velocity snapshots with shape `(nCells, 3, nTimes)`. The `cellVolumes` files are used to weight the POD inner product.

---

## 3. Compute POD Modes and Low-Order Reconstruction

The helper script `tutorials/compute_pod.py` reads the per-processor `foamToNumpy` output, concatenates pressure plus the first two velocity components into a global state matrix, weights each cell by `sqrt(cellVolume)`, computes the singular value decomposition, removes the weighting from the reconstructed fields and modes, and splits the results back into the original processor partitions.

Run:

```bash
cd tutorials
python3 compute_pod.py \
    --input of_cavity/exported_data \
    --reconstruction-output reconstruction_data \
    --modes-output mode_data \
    --rank 15 \
    --n-modes 5
```

Internally, the POD computation is:

```python
state = np.vstack((p_snap, U_snap[:, 0, :], U_snap[:, 1, :]))
weights = np.sqrt(np.tile(cell_volumes, 3))[:, None]
weighted_state = weights * state

U, D, VT = np.linalg.svd(weighted_state, full_matrices=False)

rank = 15
weighted_reconstruction = U[:, :rank] @ np.diag(D[:rank]) @ VT[:rank, :]
state_reconstructed = weighted_reconstruction / weights
modes = U[:, :5] / weights
```

The script writes two NumPy datasets in the layout expected by `numpyToFoam`:

```text
reconstruction_data/
├── p/
│   ├── p_proc_0.npy
│   ├── p_proc_1.npy
│   └── ...
└── U/
    ├── U_proc_0.npy
    ├── U_proc_1.npy
    └── ...

mode_data/
├── p_modes/
│   ├── p_modes_proc_0.npy
│   ├── p_modes_proc_1.npy
│   └── ...
└── U_modes/
    ├── U_modes_proc_0.npy
    ├── U_modes_proc_1.npy
    └── ...
```

For reconstructed fields, each column is one physical time snapshot. For POD modes, each column is one mode and is later written as a separate OpenFOAM time directory. The velocity mode files contain the first two velocity components from the POD state; the third component is written as zero.

---

## 4. Write the Reconstructed Fields Back to OpenFOAM

Create a copy of the original case and clean the processor data so that only the mesh and `0/` fields remain:

```bash
cd tutorials
cp -r of_cavity of_cavity_reconstructed
cd of_cavity_reconstructed
./Clean_proc_data
cp -r ../reconstruction_data .
```

Use the following `system/numpyToFoamDict`:

```plaintext
dataDir       reconstruction_data;
fields        (p U);

time
{
    startTime   0.01;
    endTime     0.5;
    deltaT      0.01;
}
```

Now run:

```bash
NPROC="$(find . -maxdepth 1 -type d -name 'processor*' | wc -l)"
mpirun -np "$NPROC" numpyToFoam -parallel
```

This writes the reconstructed pressure and velocity snapshots into the corresponding OpenFOAM time directories.

Create a ParaView marker file:

```bash
touch of_cavity_reconstructed.foam
```

---

## 5. Write the POD Modes to OpenFOAM

The POD modes can be written in the same way by creating a separate case:

```bash
cd tutorials
cp -r of_cavity of_cavity_modes
cd of_cavity_modes
./Clean_proc_data
cp -r ../mode_data .
```

Use the following `system/numpyToFoamDict`:

```plaintext
dataDir       mode_data;
fields        (p_modes U_modes);

time
{
    startTime   1;
    endTime     5;
    deltaT      1;
}
```

Since POD modes do not have a physical time dimension, the time settings are used here only as labels. Each written time directory corresponds to one POD mode:

- time `1` -> mode 1
- time `2` -> mode 2
- time `3` -> mode 3
- and so on

Run:

```bash
NPROC="$(find . -maxdepth 1 -type d -name 'processor*' | wc -l)"
mpirun -np "$NPROC" numpyToFoam -parallel
```

The POD modes are now available as OpenFOAM fields and can be visualized directly in ParaView.

Create a ParaView marker file:

```bash
touch of_cavity_modes.foam
```

---

## Summary

This tutorial shows how to:

- export OpenFOAM snapshots to `.npy` files with `foamToNumpy`
- compute low-order reconstructed fields from pressure and in-plane velocity components
- write reconstructed fields and POD modes back into OpenFOAM format with `numpyToFoam`
