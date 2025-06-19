struct PMGSolver{T}
    ml::MultiLevel
    b::Vector{T}
end


struct SmoothedAggregationCoarseSolver{TK,TKW} <: CoarseSolver
    args::TK
    kwargs::TKW
end

struct RugeStubenCoarseSolver{TK,TKW} <: CoarseSolver 
    args::TK
    kwargs::TKW
end

struct AMGCoarseSolver{TA,TG<:AMGAlg,TK,TKW} <: CoarseSolver
    A::TA
    alg::TG
    args::TK
    kwargs::TKW
end

function (sa::SmoothedAggregationCoarseSolver)(A)
    return AMGCoarseSolver(A, SmoothedAggregationAMG(), sa.args...; sa.kwargs...)
end

function (rs::RugeStubenCoarseSolver)(A)
    return AMGCoarseSolver(A, RugeStubenAMG(), rs.args...; rs.kwargs...)
end

function (amg::AMGCoarseSolver)(x::Vector, b::Vector)
    x_amg = AMG.solve(amg.A, b, amg.alg, amg.args...; amg.kwargs...)
    x .= x_amg
end

function AMGCoarseSolver(A, alg::AMGAlg, args...; kwargs...)
    return AMGCoarseSolver(A, alg, args, kwargs)
end


function solve(A::AbstractMatrix, b::Vector, fe_space::FESpace, coarse_solver::Type{<:CoarseSolver} = AMGCoarseSolver, args...; kwargs...)
    solver = init(A, b, fe_space, coarse_solver, args...; kwargs...)
    solve!(solver, args...; kwargs...)
end

function init(A, b, fine_fespace::FESpace , pcoarse_solver = SmoothedAggregationCoarseSolver, args...; kwargs...)
    PMGSolver(pmultigrid(A, fine_fespace,coarse_solver = _setup_caorse_solver(pcoarse_solver,args...;kwargs...)), b)
end

function solve!(solt::PMGSolver, args...; kwargs...)
    _solve(solt.ml, solt.b, args...; kwargs...)
end

_setup_caorse_solver(solvertype ,args...;kwargs...) = solvertype
_setup_caorse_solver(solvertype::Type{<:SmoothedAggregationCoarseSolver}, args...; kwargs...) =  SmoothedAggregationCoarseSolver(args..., kwargs...)
_setup_caorse_solver(solvertype::Type{<:RugeStubenCoarseSolver}, args...; kwargs...) =  RugeStubenCoarseSolver(args..., kwargs...)
