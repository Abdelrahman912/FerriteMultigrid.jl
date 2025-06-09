module PolynomialMultigrid

using Ferrite
import Ferrite: getorder
using AlgebraicMultigrid
import AlgebraicMultigrid:
    AMGAlg, Level, CoarseSolver, Multilevel, residual!, coarse_x!, coarse_b!, Pinv


include("fe.jl")
include("prolongator.jl")
include("pmultigrid.jl")
include("multilevel.jl")

end
