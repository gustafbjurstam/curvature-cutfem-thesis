# Code for *Cut Finite Element Approximation of the Mean Curvature Vector*

This repository contains the code used to generate the numerical data and figures in my master's thesis *Cut Finite Element Approximation of the Mean Curvature Vector*. You can find more details on the [projects tab of my website](https://bjurstam.eu/projects).

The code depends on the upstream [`CutFEM-Library`](https://github.com/CutFEM/CutFEM-Library/tree/development). For reproducibility, the thesis computations were run with the `development` branch at commit `952105b`.

## Requirements
A working setup needs:

- a C++20-capable compiler;
- CMake 3.16 or newer;
- MPI;
- UMFPACK/SuiteSparse for the default CutFEM configuration;
- LAPACKE and LAPACK for the Algoim-based targets;
- the upstream [`CutFEM-Library`](https://github.com/CutFEM/CutFEM-Library/tree/development) checked out at commit `952105b` on the `development` branch.

Build `CutFEM-Library` using its own instructions. This project assumes that the library has already been built and that its location is passed through `CUTFEM_ROOT`.

## Building

Create a build directory and configure the project with the path to the CutFEM library:

```bash
mkdir build
cd build
cmake -DCUTFEM_ROOT=/path/to/CutFEM-Library/ ..
```

Then build one of the executable targets, for example:

```bash
make curvature_2d
```

Available targets are:

```bash
make curvature_2d
make curvature_2d_algoim
make curvature_3d
make curvature_3d_algoim
```

If the CutFEM build directory is not located at `/path/to/CutFEM-Library/build`, pass it explicitly:

```bash
cmake \
  -DCUTFEM_ROOT=/path/to/CutFEM-Library/ \
  -DCUTFEM_BUILD_DIR=/path/to/CutFEM-Library/build \
  ..
```

## Runtime and memory requirements

The full reproduction runs are computationally expensive. On my workstation, equipped with an AMD Ryzen 9 5900X and 32 GB RAM, most full parameter sweeps were typically started in the evening and finished around the following morning.

The Algoim-based runs are particularly expensive, with quadrature being the dominant cost. Although this part of the computation is highly parallel in principle, memory usage limited the number of MPI ranks I could use in practice.

For this reason, it is recommended to first build and run a single target as a smoke test before launching the full reproduction scripts.

## Smoke test

After building, a single example can be run from the build directory, for instance:

```bash
mpirun -x OMP_NUM_THREADS=8 -n 1 ./bin/curvature_2d
```
This checks that the build and runtime setup work, but it does not reproduce all thesis data.

## Running the experiments

The executables are configured at compile time through marked `#define` lines in the source files. The scripts in `scripts/` automate this process by editing those definitions, rebuilding the corresponding target, and running the executable.

These scripts are intended for reproducing the thesis parameter sweeps, not for quick testing. The scripts are deliberately simple and can be edited to better match the available hardware. In particular, users may want to change the number of MPI ranks, the value of `OMP_NUM_THREADS`, or restrict the loops to only a subset of geometries, meshes, or stabilization choices.

The scripts assume that they are launched from the `build` directory created above. For example:

```bash
cd build
bash ../scripts/auto_2d.sh
```
The currently provided reproduction scripts are:

```bash
bash ../scripts/auto_2d.sh
bash ../scripts/auto_3d.sh
bash ../scripts/auto_2d_algoim.sh
bash ../scripts/auto_3d_algoim.sh
```

The scripts correspond to different groups of experiments:
* `auto_2d.sh` runs the two-dimensional $\mathbb{P}_1$/$\mathbb{Q}_1$ experiments;
* `auto_3d.sh` runs the three-dimensional $\mathbb{P}_1$ experiments;
* `auto_2d_algoim.sh` runs the two-dimensional Algoim experiments;
* `auto_3d_algoim.sh` runs the three-dimensional Algoim experiments.

The scripts modify marked `#define` lines in the source files. After running them, the working tree may therefore contain modified source files. To reset these changes, run from the repository root:

```bash
git restore src/curvature/*.cpp
```

## Post-processing and plotting

The MATLAB scripts in `matlab/` are used to process the generated data and produce the convergence and condition-number plots. They are intended to be run after the corresponding C++ experiments have produced output under `results/`.

Run the MATLAB scripts from the `matlab/` directory.

The main plotting scripts are:

```bash
matlab/convergence_2d.m
matlab/convergence_3d.m
matlab/convergence_2d_algoim.m
matlab/convergence_3d_algoim.m
matlab/cond_2d.m
matlab/cond_2d_algoim.m
```

The convergence scripts call `process_example_data.m`, while the condition-number scripts call `process_condition_number.m`. The lower-level plotting is handled by helper functions in the same directory.

The scripts assume that the generated data are stored under `../results/`, relative to the `matlab/` directory. For example, the expected output directories include:

```text
results/curvature_p1_2d
results/curvature_p1_3d
results/curvature_2d_algoim
results/curvature_3d_algoim
```

The MATLAB scripts are deliberately short and can be edited directly to select which examples to process, change plot ranges, or adapt paths to a different local directory layout.
