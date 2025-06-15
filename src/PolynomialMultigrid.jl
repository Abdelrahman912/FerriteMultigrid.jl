module PolynomialMultigrid

using Reexport
import Base: *

using Ferrite
import Ferrite: getorder, AbstractDofHandler, AbstractCell, AbstractRefShape
@reexport using AlgebraicMultigrid
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
include("multigrid_problems.jl")
include("prolongator.jl")
include("pmultigrid.jl")
include("multilevel.jl")
include("gallery.jl")

export FESpace, DiffusionMultigrid, ConstantCoefficient

end
