#!/bin/sh
set -e

cd "${0%/*}" || exit 1

(cd ../.. && ./Allwmake)
./reset_tutorial.sh

echo "Running the cavity with online NumPy export"
(cd of_cavity && ./Allrun)

echo "Computing POD data"
python3 compute_pod.py

echo "Importing reconstructed snapshots"
cp -a of_cavity of_cavity_reconstructed
(cd of_cavity_reconstructed && ./Clean_proc_data)
cp -a reconstruction_data of_cavity_reconstructed/
(cd of_cavity_reconstructed && \
    mpirun -np 4 numpyPostProcess -parallel \
        -dict numpyPostProcessReconstructionDict)
touch of_cavity_reconstructed/of_cavity_reconstructed.foam

echo "Importing POD modes"
cp -a of_cavity of_cavity_modes
(cd of_cavity_modes && ./Clean_proc_data)
cp -a mode_data of_cavity_modes/
(cd of_cavity_modes && \
    mpirun -np 4 numpyPostProcess -parallel \
        -dict numpyPostProcessModesDict)
touch of_cavity_modes/of_cavity_modes.foam

echo "Tutorial complete"
