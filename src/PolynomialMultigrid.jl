module PolynomialMultigrid

using Ferrite
import Ferrite: getorder, AbstractDofHandler
using AlgebraicMultigrid
import AlgebraicMultigrid:
    AMGAlg,
    Level,
    CoarseSolver,
    MultiLevel,
    residual!,
    coarse_x!,
    coarse_b!,
    Pinv,
    MultiLevelWorkspace
using LinearAlgebra
using SparseArrays

include("fe.jl")
include("prolongator.jl")
include("pmultigrid.jl")
include("multilevel.jl")

export FESpace

end
