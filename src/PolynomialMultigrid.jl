module PolynomialMultigrid

using Ferrite
using AlgebraicMultigrid
import AlgebraicMultigrid: AMGAlg, Level, CoarseSolver, Multilevel

include("prolongator.jl")
include("pmultigrid.jl")
include("multilevel.jl")

end
