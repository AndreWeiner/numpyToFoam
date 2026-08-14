#!/bin/bash
#------------------------------------------------------------------------------

# It pulls prebuilt Apptainer images from GHCR and runs the numpyToFoam tests.

# definition of path variables
main_folder="$GITHUB_WORKSPACE/unittest"
echo "$main_folder"

run_folder="$main_folder/run"
cavity_prepare="$main_folder/prepare_cavity"
versions_folder_name="$main_folder/of_versions"
SRC_FOLDER="$GITHUB_WORKSPACE/src"

echo "Starting GitHub Actions test"

mkdir -p "$run_folder"
mkdir -p "$versions_folder_name"

overall_ok=1

# ----------------------------------------------------------------------
# function: make checksum manifest for processor files
# ----------------------------------------------------------------------
makeChecksumManifest() {
    case_dir="$1"
    out_file="$2"

    (
        cd "$case_dir" || exit 1

        find processor* -type f \
            ! -path "*/constant/polyMesh/*" \
            ! -path "*/uniform/*" \
            ! -name "phi" \
            -print0 \
        | sort -z \
        | xargs -0 -r md5sum
    ) > "$out_file"
}

# ----------------------------------------------------------------------
# function: compare two checksum manifests
# ----------------------------------------------------------------------
compareChecksumManifest() {
    ref_file="$1"
    new_file="$2"
    diff_file="$3"

    ref_hashes="${ref_file}.hashes"
    new_hashes="${new_file}.hashes"

    awk '{print $1}' "$ref_file" > "$ref_hashes"
    awk '{print $1}' "$new_file" > "$new_hashes"

    diff -u "$ref_hashes" "$new_hashes" > "$diff_file" 2>&1
    status=$?

    rm -f "$ref_hashes" "$new_hashes"
    return $status
}

# ----------------------------------------------------------------------
# OpenFOAM versions to test
# ----------------------------------------------------------------------
versions=(2506 2512)

# ----------------------------------------------------------------------
# Pull Apptainer images from GitHub Container Registry
# ----------------------------------------------------------------------
for version in "${versions[@]}"; do
    sif_file="$versions_folder_name/openfoam$version.sif"

    if [ ! -f "$sif_file" ]; then
        echo "Pulling $sif_file from GHCR"

        apptainer pull "$sif_file" \
            "oras://ghcr.io/tanujravi/openfoam$version:default"

        if [ "$?" -ne 0 ]; then
            echo "Failed to pull openfoam$version.sif"
            overall_ok=0
        fi
    else
        echo "$sif_file already exists, skipping pull"
    fi
done

# ----------------------------------------------------------------------
# Run tests for each OpenFOAM version
# ----------------------------------------------------------------------
for version in "${versions[@]}"; do

    foamToNumpy_build_ok=0
    function_object_build_ok=0
    numpy_post_process_build_ok=0
    numpyToFoam_build_ok=0
    allrun_ok=0
    function_object_run_ok=0
    function_object_check_ok=0
    foamToNumpy_run_ok=0
    clean_proc_data_ok=0
    numpyToFoam_run_ok=0
    checksum_ok=0
    function_object_suite_ok=0

    cd "$run_folder" || exit 1

    version_folder="$run_folder/of$version"
    rm -rf "$version_folder"
    mkdir -p "$version_folder"
    cd "$version_folder" || exit 1

    image="$versions_folder_name/openfoam$version.sif"
    source="source /usr/lib/openfoam/openfoam$version/etc/bashrc"

    log_foamToNumpy_build="$version_folder/log.foamToNumpy.build"
    log_function_object_build="$version_folder/log.functionObject.build"
    log_numpy_post_process_build="$version_folder/log.numpyPostProcess.build"
    log_numpyToFoam_build="$version_folder/log.numpyToFoam.build"
    log_allrun="$version_folder/log.Allrun"
    log_function_object_run="$version_folder/log.functionObject"
    log_foamToNumpy_run="$version_folder/log.foamToNumpy"
    log_clean_proc_data="$version_folder/log.Clean_proc_data"
    log_numpyToFoam_run="$version_folder/log.numpyToFoam"
    log_function_object_suite="$version_folder/log.functionObjectSuite"

    checksum_ref="$version_folder/checksum.reference.md5"
    checksum_final="$version_folder/checksum.final.md5"
    checksum_diff="$version_folder/checksum.diff"

    : > "$log_foamToNumpy_build"
    : > "$log_function_object_build"
    : > "$log_numpy_post_process_build"
    : > "$log_numpyToFoam_build"
    : > "$log_allrun"
    : > "$log_function_object_run"
    : > "$log_foamToNumpy_run"
    : > "$log_clean_proc_data"
    : > "$log_numpyToFoam_run"
    : > "$log_function_object_suite"

    if [ ! -f "$image" ]; then
        echo "Image not found: $image"
        overall_ok=0
        continue
    fi

    echo "----------------------------------------"
    echo "Testing OpenFOAM version: $version"
    echo "Using image: $image"
    echo "----------------------------------------"

    # Check OpenFOAM environment inside the image
    apptainer exec "$image" bash -lc "$source && echo OpenFOAM version: \$WM_PROJECT_VERSION && command -v blockMesh"
    if [ "$?" -ne 0 ]; then
        echo "OpenFOAM environment check failed for version $version"
        overall_ok=0
        continue
    fi

    cp -r "$SRC_FOLDER" .

    FOAMTONUMPY="$version_folder/src/foamToNumpy"
    cd "$FOAMTONUMPY" || exit 1
    apptainer exec "$image" bash -lc "$source && wmake" \
        > "$log_foamToNumpy_build" 2>&1
    foamToNumpy_build_ok=$(( $? == 0 ))

    FUNCTION_OBJECT="$version_folder/src/functionObjects/foamToNumpy"
    cd "$FUNCTION_OBJECT" || exit 1
    apptainer exec "$image" bash -lc "$source && wmake" \
        > "$log_function_object_build" 2>&1
    function_object_build_ok=$(( $? == 0 ))

    NUMPY_POST_PROCESS="$version_folder/src/numpyPostProcess"
    cd "$NUMPY_POST_PROCESS" || exit 1
    apptainer exec "$image" bash -lc "$source && wmake" \
        > "$log_numpy_post_process_build" 2>&1
    numpy_post_process_build_ok=$(( $? == 0 ))

    NUMPYTOFOAM="$version_folder/src/numpyToFoam"
    cd "$NUMPYTOFOAM" || exit 1
    apptainer exec "$image" bash -lc "$source && wmake" \
        > "$log_numpyToFoam_build" 2>&1
    numpyToFoam_build_ok=$(( $? == 0 ))

    case_dir="$version_folder/of_cavity"
    log_icoFoam="$case_dir/log.icoFoam"
    apptainer exec "$image" bash -lc \
        "$source && '$cavity_prepare' '$case_dir' && cd '$case_dir' && export OMPI_MCA_rmaps_base_oversubscribe=1 && ./Allrun" \
        > "$log_allrun" 2>&1
    allrun_ok=$(( $? == 0 ))
    cd "$case_dir" || exit 1

    # save reference checksums after original simulation
    if [ "$allrun_ok" -eq 1 ]; then
        makeChecksumManifest "$case_dir" "$checksum_ref"
    fi

    apptainer exec "$image" bash -lc \
        "$source && mpirun --oversubscribe -np 4 postProcess -parallel -dict system/controlDict.numpy -time '0.1:0.5'" \
        > "$log_function_object_run" 2>&1
    function_object_run_ok=$(( $? == 0 ))

    if [ "$function_object_run_ok" -eq 1 ]; then
        python3 "$main_folder/check_function_object.py" \
            "$case_dir/postProcessing/numpyExport" 0.1 0.2 0.3 0.4 0.5
        function_object_check_ok=$(( $? == 0 ))
    fi

    apptainer exec "$image" bash -lc "$source && mpirun --oversubscribe -np 4 foamToNumpy -parallel" \
        > "$log_foamToNumpy_run" 2>&1
    foamToNumpy_run_ok=$(( $? == 0 ))

    apptainer exec "$image" bash -lc "$source && ./Clean_proc_data" \
        > "$log_clean_proc_data" 2>&1
    clean_proc_data_ok=$(( $? == 0 ))

    apptainer exec "$image" bash -lc "$source && mpirun --oversubscribe -np 4 numpyToFoam -parallel" \
        > "$log_numpyToFoam_run" 2>&1
    numpyToFoam_run_ok=$(( $? == 0 ))

    # compare checksums after numpyToFoam
    if [ "$numpyToFoam_run_ok" -eq 1 ] && [ -f "$checksum_ref" ]; then
        makeChecksumManifest "$case_dir" "$checksum_final"

        if compareChecksumManifest "$checksum_ref" "$checksum_final" "$checksum_diff"; then
            checksum_ok=1
        else
            checksum_ok=0
        fi
    fi

    apptainer exec "$image" bash -lc \
        "$source && export OMPI_MCA_rmaps_base_oversubscribe=1 && cd '$main_folder' && ./Allrun-functionObjects" \
        > "$log_function_object_suite" 2>&1
    function_object_suite_ok=$(( $? == 0 ))
    
    foamToNumpy_build_status=$([ "$foamToNumpy_build_ok" -eq 1 ] && echo "success" || echo "failed")
    function_object_build_status=$([ "$function_object_build_ok" -eq 1 ] && echo "success" || echo "failed")
    numpy_post_process_build_status=$([ "$numpy_post_process_build_ok" -eq 1 ] && echo "success" || echo "failed")
    numpyToFoam_build_status=$([ "$numpyToFoam_build_ok" -eq 1 ] && echo "success" || echo "failed")
    allrun_status=$([ "$allrun_ok" -eq 1 ] && echo "success" || echo "failed")
    function_object_run_status=$([ "$function_object_run_ok" -eq 1 ] && echo "success" || echo "failed")
    function_object_check_status=$([ "$function_object_check_ok" -eq 1 ] && echo "success" || echo "failed")
    foamToNumpy_run_status=$([ "$foamToNumpy_run_ok" -eq 1 ] && echo "success" || echo "failed")
    clean_proc_data_status=$([ "$clean_proc_data_ok" -eq 1 ] && echo "success" || echo "failed")
    numpyToFoam_run_status=$([ "$numpyToFoam_run_ok" -eq 1 ] && echo "success" || echo "failed")
    checksum_status=$([ "$checksum_ok" -eq 1 ] && echo "success" || echo "failed")
    function_object_suite_status=$([ "$function_object_suite_ok" -eq 1 ] && echo "success" || echo "failed")
    
    echo "----------------------------------------"
    echo "Version: $version"
    echo "foamToNumpy build      : $foamToNumpy_build_status"
    echo "function object build  : $function_object_build_status"
    echo "numpyPostProcess build : $numpy_post_process_build_status"
    echo "numpyToFoam build      : $numpyToFoam_build_status"
    echo "Allrun                 : $allrun_status"
    echo "function object run    : $function_object_run_status"
    echo "function object check  : $function_object_check_status"
    echo "foamToNumpy run        : $foamToNumpy_run_status"
    echo "Clean_proc_data        : $clean_proc_data_status"
    echo "numpyToFoam run        : $numpyToFoam_run_status"
    echo "checksum match         : $checksum_status"
    echo "full function suite    : $function_object_suite_status"
    echo "----------------------------------------"
    
    if [ "$foamToNumpy_build_ok" -ne 1 ] || \
       [ "$function_object_build_ok" -ne 1 ] || \
       [ "$numpy_post_process_build_ok" -ne 1 ] || \
       [ "$numpyToFoam_build_ok" -ne 1 ] || \
       [ "$allrun_ok" -ne 1 ] || \
       [ "$function_object_run_ok" -ne 1 ] || \
       [ "$function_object_check_ok" -ne 1 ] || \
       [ "$foamToNumpy_run_ok" -ne 1 ] || \
       [ "$clean_proc_data_ok" -ne 1 ] || \
       [ "$numpyToFoam_run_ok" -ne 1 ] || \
       [ "$checksum_ok" -ne 1 ] || \
       [ "$function_object_suite_ok" -ne 1 ]; then

        echo "Failure detected for OpenFOAM version $version"
        echo "Showing useful logs"

        echo "----- log.foamToNumpy.build -----"
        tail -n 80 "$log_foamToNumpy_build" || true

        echo "----- log.functionObject.build -----"
        tail -n 80 "$log_function_object_build" || true

        echo "----- log.numpyPostProcess.build -----"
        tail -n 80 "$log_numpy_post_process_build" || true

        echo "----- log.numpyToFoam.build -----"
        tail -n 80 "$log_numpyToFoam_build" || true

        echo "----- log.Allrun -----"
        tail -n 80 "$log_allrun" || true

        echo "----- Allrun file -----"
        cat "$case_dir/Allrun" || true
        
        echo "----- log.icoFoam -----"
        tail -n 80 "$log_icoFoam" || true

        echo "----- log.functionObject -----"
        tail -n 80 "$log_function_object_run" || true
        
        echo "----- log.foamToNumpy -----"
        tail -n 80 "$log_foamToNumpy_run" || true

        echo "----- log.Clean_proc_data -----"
        tail -n 80 "$log_clean_proc_data" || true

        echo "----- log.numpyToFoam -----"
        tail -n 80 "$log_numpyToFoam_run" || true

        echo "----- checksum.diff -----"
        cat "$checksum_diff" || true

        echo "----- log.functionObjectSuite -----"
        tail -n 160 "$log_function_object_suite" || true

        overall_ok=0
    fi

done

if [ "$overall_ok" -eq 1 ]; then
    echo "All tests passed"
    exit 0
else
    echo "One or more tests failed"
    exit 1
fi
