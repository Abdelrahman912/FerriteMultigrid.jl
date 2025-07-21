# FerriteMultigrid.jl

[![Build Status](https://github.com/abdelrahman912/FerriteMultigrid.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/abdelrahman912/FerriteMultigrid.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/abdelrahman912/FerriteMultigrid.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/abdelrahman912/FerriteMultigrid.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://abdelrahman912.github.io/FerriteMultigrid.jl/dev/)


**FerriteMultigrid.jl** is a lightweight, flexible **p-multigrid framework** designed for high-order finite element problems in Julia.  
It is built on top of [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) and leverages [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) as the coarse-grid solver once the approximation is reduced to \( p = 1 \).


## Example Usage

```julia
using FerriteMultigrid

# Define a 1D diffusion problem with p = 2 and 3 quadrature points.
K, f, fe_space = poisson(1000, 2, 3)

# Define a p-multigrid configuration
config = pmultigrid_config() # default config (galerkin as coarsening strategy and direct projection (i.e., from p to 1 directly))

# Solve using the p-multigrid solver
x, res = solve(K, f, fe_space, config; log = true, rtol = 1e-10)
```

## Acknowledgement

This framework is primarily developed at the [chair of continuum mechanics at Ruhr University Bochum](https://www.lkm.ruhr-uni-bochum.de/) under 
the supervision of @termi-official.
