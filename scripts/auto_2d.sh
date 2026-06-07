#!/bin/bash

# Path to the C++ file (update this if you are running the script from a build folder)
CPP_FILE="../src/curvature/curvature_2d.cpp"

GEOMETRIES=("PEANUT" "STAR" "AMOEBA")
# P1
STABILISATIONS=("elem" "patch")

sed -i "s/^#define .* \/\/ SCRIPT_MARKER2/#define TRIANGLES \/\/ SCRIPT_MARKER2/" "$CPP_FILE"

for geo in "${GEOMETRIES[@]}"; do
    for stab in "${STABILISATIONS[@]}"; do
        TARGET="${geo}_${stab}"
        echo "========================================"
        echo "Configuring and compiling for: $TARGET"
        echo "========================================"
        
        # Hot-swap the #define line using the SCRIPT_MARKER
        # Note: If you are on macOS, you need to use `sed -i '' "s/^#define.*/..."` instead
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER1/#define ${TARGET} \/\/ SCRIPT_MARKER1/" "$CPP_FILE"
        
        # Compile (adjust this to your specific make target if necessary)
        make curvature_2d
        
        # Run the executable
        # Adjust the MPI cores or executable path as needed for your environment
        mpirun -x OMP_NUM_THREADS=8 -n 1 ./bin/curvature_2d
        
    done
done

# Q1
sed -i "s/^#define .* \/\/ SCRIPT_MARKER2/#define QUADS \/\/ SCRIPT_MARKER2/" "$CPP_FILE"
STABILISATIONS=("elem" "surf")

for geo in "${GEOMETRIES[@]}"; do
    for stab in "${STABILISATIONS[@]}"; do
        TARGET="${geo}_${stab}"
        echo "========================================"
        echo "Configuring and compiling for: $TARGET"
        echo "========================================"
        
        # Hot-swap the #define line using the SCRIPT_MARKER
        # Note: If you are on macOS, you need to use `sed -i '' "s/^#define.*/..."` instead
        sed -i "s/^#define .* \/\/ SCRIPT_MARKER1/#define ${TARGET} \/\/ SCRIPT_MARKER1/" "$CPP_FILE"
        
        # Compile (adjust this to your specific make target if necessary)
        make curvature_2d
        
        # Run the executable
        # Adjust the MPI cores or executable path as needed for your environment
        mpirun -x OMP_NUM_THREADS=8 -n 1 ./bin/curvature_2d
        
    done
done



echo "All configurations completed."