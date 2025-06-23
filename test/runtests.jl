using FerriteMultigrid
import FerriteMultigrid: _element_mass_matrix!, _element_prolongator!, build_prolongator
using Test
using Ferrite
using SparseArrays


include("test_prolongator.jl")
include("test_multilevel.jl")
include("test_examples.jl")
