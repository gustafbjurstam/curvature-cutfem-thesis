#!/bin/bash

# Path to the C++ file (update this if you are running the script from a build folder)
CPP_FILE="../src/curvature/curvature_3d_algoim.cpp"

MESHES=("TET_P2" "HEXA_Q1")
GEOMETRIES=("ELL_TOR" "DECO_CUBE")

for mesh in "${MESHES[@]}"; do
    for geo in "${GEOMETRIES[@]}"; do
        TARGET="${mesh}_${geo}_ELEM_STAB"
        echo "========================================"
        echo "Configuring and compiling for: $TARGET"
        echo "========================================"

        # Hot-swap the #define lines using the SCRIPT_MARKERs
        # Note: If you are on macOS, you need to use `sed -i '' "s/^#define.*/.../"`
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER1/#define ${geo} \/\/ SCRIPT_MARKER1/" "$CPP_FILE"
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER2/#define ELEM_STAB \/\/ SCRIPT_MARKER2/" "$CPP_FILE"
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER3/#define ${mesh} \/\/ SCRIPT_MARKER3/" "$CPP_FILE"

        # Compile
        make curvature_3d_algoim

        # Run the executable
        # Adjust the MPI cores or executable path as needed for your environment
        mpirun -x OMP_NUM_THREADS=12 -n 1 ./bin/curvature_3d_algoim

    done
done

echo "All configurations completed."