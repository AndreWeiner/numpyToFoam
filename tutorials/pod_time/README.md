# Time-resolved POD

This tutorial uses the NumPy function objects for a complete reduced-order
workflow:

1. `foamToNumpy` exports `p`, `U`, times, and cell volumes while `icoFoam`
   runs.
2. [compute_pod.py](compute_pod.py) computes a volume-weighted POD, a rank-15
   reconstruction, and the first five modes.
3. `numpyPostProcess` imports the generated batches directly with
   `numpyToFoam`. The importer writes native fields here because the two output
   cases are intended for visualization.

Run the tutorial from an initialized OpenFOAM environment:

```sh
cd tutorials/pod_time
./run_tutorial.sh
```

The result is two decomposed cases:

- `of_cavity_reconstructed`, containing reconstructed snapshots at the
  original physical times;
- `of_cavity_modes`, using times `1` through `5` as mode labels.

The online exporter is configured in
`of_cavity/system/controlDict`. Import is configured by the two
`numpyPostProcess` dictionaries in the same `system` directory. No legacy
`foamToNumpy` or `numpyToFoam` application is involved.

The Python code deliberately has no command-line interface. The POD rank and
number of displayed modes are constants near the top of the script, keeping
the example focused on the array operations and batch layout.

Remove generated data with:

```sh
./reset_tutorial.sh
```
