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


function solve(
    A::AbstractMatrix,
    b::Vector,
    fe_space::FESpace#=, coarse_solver::CoarseSolver =#
)
    #solver = init(A, b, fe_space, coarse_solver)
    solver = init(A, b, fe_space)
    solve!(solver)
end

function init(A, b, fine_fespace::FESpace#=, coarse_solver::CoarseSolver=#)
    #PMGSolver(pmultigrid(A, fine_fespace, coarse_solver), b)
    PMGSolver(pmultigrid(A, fine_fespace), b)
end

function solve!(solt::PMGSolver)
    _solve(solt.ml, solt.b;log=true, reltol=1e-10 )
end
