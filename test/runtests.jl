using PolynomialMultigrid
import PolynomialMultigrid: _element_mass_matrix!, _element_prolongator!, build_prolongator
using Test
using Ferrite
using Tensors
using SparseArrays


include("test_prolongator.jl")
include("test_examples.jl")
