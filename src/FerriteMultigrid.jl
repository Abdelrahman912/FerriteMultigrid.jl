module FerriteMultigrid

using Reexport
using Ferrite
import Ferrite: getorder, AbstractDofHandler
@reexport using AlgebraicMultigrid
import AlgebraicMultigrid as AMG
import AlgebraicMultigrid:
    AMGAlg,
    Level,
    CoarseSolver,
    MultiLevel,
    residual!,
    coarse_x!,
    coarse_b!,
    Pinv,
    _solve,
    MultiLevelWorkspace
using LinearAlgebra
using SparseArrays

include("fe.jl")
include("prolongator.jl")
include("pmultigrid.jl")
include("multilevel.jl")
include("gallery.jl")

export FESpace, SmoothedAggregationCoarseSolver, RugeStubenCoarseSolver

end
