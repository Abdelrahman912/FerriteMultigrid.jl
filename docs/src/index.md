```@meta
DocTestSetup = :(using FerriteMultigrid)
```

## FerriteMultigrid

**FerriteMultigrid.jl** is a lightweight, flexible **p-multigrid framework** designed for high-order finite element problems in Julia.  
It is built on top of [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) and leverages [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) as the coarse-grid solver once the approximation is reduced to \( p = 1 \).


## How the documentation is organized

This high level view of the documentation structure will help you find what you are looking
for. The document is organized as follows:

 - [**Tutorials**](tutorials/index.md) are thoroughly documented examples which guides you
   through the process of solving finite element problems using p- and A- multigrid methods.
 - [**API Reference**](api-reference/index.md) contains the technical API reference of functions and
   methods (e.g. the documentation strings).
