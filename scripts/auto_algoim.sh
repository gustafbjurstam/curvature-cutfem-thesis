#!/bin/bash

# Path to the C++ file (update this if you are running the script from a build folder)
CPP_FILE="../src/curvature/curvature_2d_algoim.cpp"

GEOMETRIES=("PEANUT" "STAR" "AMOEBA")
STABILISATIONS=("elem" "surf")

for geo in "${GEOMETRIES[@]}"; do
    for stab in "${STABILISATIONS[@]}"; do
        TARGET="${geo}_${stab}"
        echo "========================================"
        echo "Configuring and compiling for: $TARGET"
        echo "========================================"
        
        # Hot-swap the #define line using the SCRIPT_MARKER
        # Note: If you are on macOS, you need to use `sed -i '' "s/^#define.*/..."` instead
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER/#define ${TARGET} \/\/ SCRIPT_MARKER/" "$CPP_FILE"
        
        # Compile (adjust this to your specific make target if necessary)
        make curvature_2d_algoim
        
        # Run the executable
        # Adjust the MPI cores or executable path as needed for your environment
        mpirun -x OMP_NUM_THREADS=8 -n 1 ./bin/curvature_2d_algoim
        
    done
done

echo "All configurations completed."