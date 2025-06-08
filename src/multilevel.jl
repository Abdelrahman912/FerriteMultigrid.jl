struct FineFESpace{DH<:AbstractDofHandler,CV<:AbstractCellValues}
    dh::DH
    cv::CV
end

struct PMGSolver{T}
    ml::MultiLevel
    b::Vector{T}
end

struct AMGCoarseSolver{A<:AMGAlg,TK,TKW}<:CoarseSolver
    alg::A
    args::TK
    kwargs::TKW
end

function (amg::AMGCoarseSolver)(A::AbstractMatrix, b::Vector)
    solve(A, b, amg.alg, amg.args...; amg.kwargs...)
end

function AMGCoarseSolver(alg::AMGAlg, args...; kwargs...)
    return AMGCoarseSolver(alg, args, kwargs)
end


function solve(A::AbstractMatrix, b::Vector, fe_space::FineFESpace, coarse_solver::CoarseSolver)
    # This function will use the AMG algorithm to solve the system Ax = b
    # using the multigrid structure defined in PMultigrid.
    # The implementation details will depend on the specific AMG algorithm used.

    # Placeholder for actual implementation
    throw("AMG solve not implemented yet")

end

function init(A, b, fe_space::FineFESpace, coarse_solver::CoarseSolver)
    PMGSolver(pmultigrid(A,fe_space,coarse_solver), b)
end

function solve!(solt::PMGSolver) 
    _solve(solt.ml, solt.b)   
end
