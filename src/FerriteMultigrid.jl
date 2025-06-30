module FerriteMultigrid

using Reexport
using LinearAlgebra
import LinearSolve
using SparseArrays
import SparseArrays: AbstractSparseMatrixCSC
@reexport import CommonSolve: solve, solve!, init
using Ferrite
import Ferrite: getorder, AbstractDofHandler, reinit!
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

include("fe.jl")
include("prolongator.jl")
include("pmultigrid.jl")
include("multilevel.jl")
include("gallery.jl")
include("precs.jl")

export 
    FESpace, 
    SmoothedAggregationCoarseSolver, 
    RugeStubenCoarseSolver, 
    Pinv,
    PolynomialMultigridPreconBuilder

end
