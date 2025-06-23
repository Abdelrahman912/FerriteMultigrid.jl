struct PMGSolver{T}
    ml::MultiLevel
    b::Vector{T}
end

struct SmoothedAggregationCoarseSolver{TK,TKW} <: CoarseSolver
    args::TK
    kwargs::TKW
end

SmoothedAggregationCoarseSolver(args...; kwargs...) = SmoothedAggregationCoarseSolver(args, kwargs)
    
struct RugeStubenCoarseSolver{TK,TKW} <: CoarseSolver 
    args::TK
    kwargs::TKW
end

RugeStubenCoarseSolver(args...; kwargs...) = RugeStubenCoarseSolver(args, kwargs)

struct AMGCoarseSolver{TA,TG<:AMGAlg,TK,TKW} <: CoarseSolver
    A::TA
    alg::TG
    args::TK
    kwargs::TKW
end

function (sa::SmoothedAggregationCoarseSolver)(A)
    return AMGCoarseSolver(A, SmoothedAggregationAMG(), sa.args; sa.kwargs)
end

function (rs::RugeStubenCoarseSolver)(A)
    return AMGCoarseSolver(A, RugeStubenAMG(), rs.args; rs.kwargs)
end

function (amg::AMGCoarseSolver)(x::Vector, b::Vector)
    x_amg = AMG.solve(amg.A, b, amg.alg, amg.args...; amg.kwargs...)
    x .= x_amg
end

function AMGCoarseSolver(A, alg::AMGAlg, args...; kwargs...)
    return AMGCoarseSolver(A, alg, args, kwargs)
end


function solve(A::AbstractMatrix, b::Vector, fe_space::FESpace, pcoarse_solvertype::Type{<:CoarseSolver} = SmoothedAggregationCoarseSolver, args...; kwargs...)
    solver = init(A, b, fe_space, pcoarse_solvertype, args...; kwargs...)
    solve!(solver, args...; kwargs...)
end

function init(A, b, fine_fespace::FESpace , pcoarse_solvertype = SmoothedAggregationCoarseSolver, args...; kwargs...)
    PMGSolver(pmultigrid(A, fine_fespace, _setup_coarse_solver(pcoarse_solvertype,args...;kwargs...),args...;kwargs...), b)
end

function solve!(solt::PMGSolver, args...; kwargs...)
    _solve(solt.ml, solt.b, args...; kwargs...)
end

_setup_coarse_solver(solvertype ,args...;kwargs...) = solvertype
_setup_coarse_solver(solvertype::Type{<:SmoothedAggregationCoarseSolver}, args...; kwargs...) =  SmoothedAggregationCoarseSolver(args..., kwargs...)
_setup_coarse_solver(solvertype::Type{<:RugeStubenCoarseSolver}, args...; kwargs...) =  RugeStubenCoarseSolver(args..., kwargs...)
