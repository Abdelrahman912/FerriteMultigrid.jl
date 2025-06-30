struct PolynomialMultigridPreconBuilder{Tk, CS}
    fe_space::FESpace
    pcoarse_solver::CS
    blocksize::Int
    kwargs::Tk
end

function PolynomialMultigridPreconBuilder(fe_space::FESpace, pcoarse_solvertype; blocksize = 1, kwargs...)
    return PolynomialMultigridPreconBuilder(fe_space, setup_coarse_solver(pcoarse_solvertype; kwargs...), blocksize, kwargs)
end

function (b::PolynomialMultigridPreconBuilder)(A::AbstractSparseMatrixCSC, p)
    return (aspreconditioner(pmultigrid(SparseMatrixCSC(A), b.fe_space, b.pcoarse_solver, Val{b.blocksize}; b.kwargs...)), I)
end
