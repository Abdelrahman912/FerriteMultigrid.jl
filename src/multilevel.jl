struct PMGSolver{T}
    ml::MultiLevel
    b::Vector{T}
end

## FIXME: AMGCoarseSolver is not designed properly yet
# struct AMGCoarseSolver{A<:AMGAlg,TK,TKW}<:CoarseSolver
#     alg::A
#     args::TK
#     kwargs::TKW
# end

struct AMGCoarseSolver{TA}<:CoarseSolver
    A::TA
end

function (amg::AMGCoarseSolver)(x::Vector, b::Vector)
    #solve(A, b, amg.alg, amg.args...; amg.kwargs...)
    x_amg = AlgebraicMultigrid.solve(amg.A, b, SmoothedAggregationAMG())
    x .= x_amg
end

function AMGCoarseSolver(alg::AMGAlg, args...; kwargs...)
    return AMGCoarseSolver(alg, args, kwargs)
end


function solve(A::AbstractMatrix, b::Vector, fe_space::FESpace, pgrid_config::PMultigridConfiguration = pmultigrid_config())
    solver = init(A, b, fe_space, pgrid_config)
    solve!(solver)
end

function init(A, b, fine_fespace::FESpace, pgrid_config::PMultigridConfiguration)
    PMGSolver(pmultigrid(A, fine_fespace, pgrid_config), b)
end

function solve!(solt::PMGSolver)
    _solve(solt.ml, solt.b;log=true, reltol=1e-10)
end
