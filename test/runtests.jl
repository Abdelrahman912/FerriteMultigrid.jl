using FerriteMultigrid
import FerriteMultigrid: _element_mass_matrix!, _element_prolongator!, build_prolongator
using Test
using Ferrite
using SparseArrays
using LinearSolve: solve, KrylovJL_CG, LinearProblem


include("test_prolongator.jl")
include("test_parameters.jl")
include("test_examples.jl")
include("test_precs.jl")
