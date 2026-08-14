# Parametric airfoil POD interpolation

This tutorial demonstrates an in-memory post-processing workflow around the
OpenFOAM `simpleFoam/airFoil2D` case.

For eight fixed angle-of-attack and speed samples it:

1. runs the full-order airfoil case;
2. exports converged `p`, `U`, and `nuTilda` online with `foamToNumpy`;
3. constructs a volume-weighted POD from the first six samples;
4. interpolates the POD coefficients with a quadratic response surface;
5. imports all predicted fields with `numpyPostProcess`; and
6. computes forces directly on the imported fields.

The final step does not create native intermediate fields. The `simpleFoam`
environment creates transport and turbulence models, updates imported boundary
conditions, refreshes the Spalart-Allmaras wall viscosity for every prediction,
and then executes the standard `forces` function object.

Run from an initialized OpenFOAM environment:

```sh
cd tutorials/pod_parametric
./run_tutorial.sh
```

The principal output is `results/comparison.csv`, containing the full-order
and POD-interpolated drag, lift, and pitching-moment coefficients. The last two
rows are test samples; the others form the interpolation set.

Configuration is intentionally kept in OpenFOAM dictionaries:

- `of_airfoil/system/controlDict` contains the online exporter;
- `of_airfoil/system/numpyPostProcessDict` selects the `simpleFoam`
  environment and the imported-field force calculation.

The Python script contains only the fixed sampling, POD/interpolation
calculation, and batch assembly. It exposes no command-line options and depends
only on NumPy.

Remove generated cases and results with:

```sh
./reset_tutorial.sh
```
