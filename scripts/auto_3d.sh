#!/bin/bash

# Path to the C++ file (update this if you are running the script from a build folder)
CPP_FILE="../src/curvature/curvature_3d.cpp"

GEOMETRIES=("DECO_CUBE" "ELL_TOR")
STABILISATIONS=("VOL_GHOST" "FACE_GHOST")

for geo in "${GEOMETRIES[@]}"; do
    for stab in "${STABILISATIONS[@]}"; do
        echo "========================================"
        echo "Configuring and compiling for: $TARGET"
        echo "========================================"
        
        # Hot-swap the #define line using the SCRIPT_MARKER
        # Note: If you are on macOS, you need to use `sed -i '' "s/^#define.*/..."` instead
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER1/#define ${geo} \/\/ SCRIPT_MARKER1/" "$CPP_FILE"
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER2/#define ${stab} \/\/ SCRIPT_MARKER2/" "$CPP_FILE"
        
        # Compile (adjust this to your specific make target if necessary)
        make curvature_3d
        
        # Run the executable
        # Adjust the MPI cores or executable path as needed for your environment
        mpirun -x OMP_NUM_THREADS=12 -n 1 ./bin/curvature_3d
        
    done
done

echo "All configurations completed."