function pmultigrid(
    A::TA,
    fe_space::FESpace,
    presmoother = GaussSeidel(),
    postsmoother = GaussSeidel(),
    coarse_solver::CoarseSolver = Pinv, # TODO: remove Pinv and add AMGCoarseSolver
)

    levels = Vector{Level{TA,TA,Adjoint{T,TA}}}()
    w = MultiLevelWorkspace(Val{bs}, eltype(A))
    residual!(w, size(A, 1))
    fine_p = fe_space |> order

    fespaces = Vector{FESpace}()
    push!(fespaces, fe_space)

    for p = (fine_p-1):-1:1
        fine_fespace = fespaces[end]
        coarse_fespace = coarsen_order(fine_fespace, p)
        push!(fespaces, coarse_fespace)

        A = _extend_hierarchy!(
            levels,
            fine_fespace,
            coarse_fespace
        )

        coarse_x!(w, size(A, 1))
        coarse_b!(w, size(A, 1))
        residual!(w, size(A, 1))
    end
    return MultiLevel(levels, A, coarse_solver(A), presmoother, postsmoother, w)
end

function extend_hierarchy!(
    levels,
    fine_fespace::FESpace,
    coarse_fespace::FESpace,
    A
)
    P = build_prolongator(fine_fespace, coarse_fespace)
    R = P' # TODO: do we need other method to compute R?
    push!(levels, Level(A, P, R))
    A = R * A * P
    dropzeros!(A)
    return A
end
