# Code for **Cut Finite Element Approximation of the Mean Curvature Vector**

This is the code used to generate the figures presented in my master's thesis **Cut Finite Element Approximation of the Mean Curvature Vector**. In order to use this a modified version of the [CutFEM-Library](https://github.com/CutFEM/CutFEM-Library/tree/development) needs to be used. My additions can be found in my [fork](https://github.com/gustafbjurstam/CutFEM-Library/tree/development). The primary changes are added capabilities for $\mathbb{Q}_1$ hexes. 

For more details, check out [my website](https://bjurstam.eu/projects).

## Running

In order to reproduce the results, an appropriate version of the CutFEM-Library is needed. Locate where it is installed, and run

```
    cmake -DCUTFEM_ROOT=/path/to/CutFEM-Library/ ..
    make example
```