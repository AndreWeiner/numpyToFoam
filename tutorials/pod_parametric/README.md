# Parametric Airfoil PODI Tutorial

This tutorial builds a small parametric reduced-order workflow around the
OpenFOAM `simpleFoam/airFoil2D` case.  The reusable case lives in
`tutorials/pod_parametric/of_airfoil`; the tutorial does not create a new
OpenFOAM case for every parameter point.  Instead, the driver updates the
freestream velocity and force directions, runs the case, exports the converged
fields with `foamToNumpy`, and stores the generated NumPy data under
`results/cases`.  The full-order fields remain decomposed; the tutorial does
not reconstruct the final solution because `foamToNumpy` reads the processor
fields directly.

The two parameters are:

- angle of attack, sampled between `-10` and `10` degrees
- freestream speed, sampled within `+/-25%` of the original airfoil speed

All parameter points are selected by a Latin-hypercube design whose bin order is
derived from a small two-dimensional Sobol sequence.  By default, one quarter of
the cases are held out for testing and the remaining cases are used to train the
POD with interpolation model.  To avoid extrapolation, the sampled parameter
points on the convex hull are always kept in the training set; only interior
sample points may become test cases.

## Run

From the repository root:

```sh
cd tutorials/pod_parametric
./run_tutorial.sh
```

Useful controls:

```sh
N_CASES=12 POD_SPLINE_SMOOTHING=0 POD_MODES=0 ./run_tutorial.sh
```

- `N_CASES` must be at least `7`; the default is `8`.
- `POD_SPLINE_SMOOTHING` is the thin-plate spline smoothing value for
  interpolating POD and force coefficients; the default is `0`, which gives an
  exact interpolant through the training samples.
- `POD_MODES=0` keeps all modes available from the training set.  Set it to a
  positive integer to truncate the basis.

The Python dependencies are NumPy for the reduced-order model and Matplotlib for
the pressure-distribution plots.

To remove generated data:

```sh
./reset_tutorial.sh
```

## What The Script Does

For every parameter point, `parametric_airfoil.py`:

1. writes the freestream vector to `of_airfoil/0.orig/U`
2. updates `magUInf`, `dragDir`, and `liftDir` in
   `of_airfoil/system/controlDict`
3. runs `of_airfoil/Allclean` and `of_airfoil/Allrun`
4. removes any previous `of_airfoil/exported_data` directory
5. exports converged `p` and `U` volume fields with `foamToNumpy`
6. archives force coefficients and sampled airfoil surface pressure

The POD state is

```text
[p, Ux, Uy]^T
```

assembled from the per-processor NumPy files written by `foamToNumpy`.  Cell
volumes are used as POD weights.  Before the POD, pressure is nondimensionalized
with the reference dynamic pressure `U_ref^2`, and velocity with `U_ref`.  Each
physical block is then min-max scaled using only the training snapshots.  PODI
predictions are inverse-scaled before `numpyToFoam` writes them back to
OpenFOAM.  The POD and force coefficients are interpolated from normalized angle
and speed with thin-plate splines, a scattered-data spline method that works
directly on the Latin/Sobol training points.

For each parameter case, including both training and test cases, the script
predicts the POD coefficients, reconstructs `p` and `U`, writes the predicted
fields back to an OpenFOAM case with `numpyToFoam`, samples the predicted airfoil
surface pressure, and compares the prediction with the full-order simulation in
normalized-x pressure-coefficient plots.  Predicted fields are sampled in
parallel as well, without a reconstructed single-processor case.  The PODI
pressure sampling uses
`simpleFoam -postProcess` so the pressure field is loaded in the solver context.
The PODI case copy preserves the `system/` dictionaries and original field
templates, but excludes solved time folders, flux files, logs, and generated
post-processing/export data.  After `numpyToFoam` writes the PODI fields, the
decomposed time-`0` `nut` and `nuTilda` fields are copied into the synthetic PODI
times so `simpleFoam -postProcess` can construct the turbulence model.  These
are auxiliary post-processing fields here; they are not part of the POD state.
Volume fields are written in OpenFOAM binary format.

## Outputs

Important generated files are:

- `results/parameters.csv`: sampled parameters and train/test labels
- `results/cases/run_*/exported_data`: `foamToNumpy` output for each full-order run
- `results/cases/run_*/surface_pressure.csv`: full-order airfoil pressure samples
- `results/cases/run_*/surface_pressure_podi.csv`: PODI airfoil pressure samples
- `results/cases/run_*/surface_pressure_coefficient_comparison.csv`: matched
  OpenFOAM/PODI `-Cp` values over normalized `x`
- `results/cases/run_*/surface_pressure_raw`: copied raw OpenFOAM pressure samples
- `results/cases/run_*/surface_pressure_podi_raw`: copied raw PODI pressure samples
- `prediction_data`: PODI-predicted NumPy fields in `numpyToFoam` layout
- `of_airfoil_podi`: OpenFOAM case containing the PODI-predicted fields
- `results/comparison.csv`: force-coefficient comparison for all cases, with a
  train/test label
- `results/scaling_reference.csv`: reference speed and pressure values
- `results/surface_pressure_plots/*.png`: OpenFOAM and PODI `-Cp` plotted over
  normalized airfoil `x`, plus their difference
- `results/force_coefficient_plots/*.png`: force-coefficient comparison plots
- `results/singular_values.npy`: POD singular values

`results/comparison.csv` contains true and predicted `Cd`, `Cl`, and pitching
moment, a train/test label, plus a relative volume-field error.  The sampled
negative pressure coefficient is compared visually over normalized `x` in one
PNG file per case and saved as a matched CSV table; points are matched by sampled
surface coordinates so parallel output order does not matter.  Since `p` is
incompressible kinematic pressure, the plot uses `Cp = p / (0.5 U_inf^2)`.

## Notes

`foamToNumpy` and `numpyToFoam` currently operate on volume fields.  This is why
the POD state is built from volume `p` and `U`, while force coefficients and
surface pressure are obtained with OpenFOAM function objects.

The implementation is intentionally small and educational: the numerical work
uses NumPy, and Matplotlib is used only to write the comparison figures.
