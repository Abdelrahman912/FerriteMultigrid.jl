function pmultigrid(
    A::TA,
    fe_space::FESpace,
    pcoarse_solver , 
    ::Type{Val{bs}} = Val{1};
    presmoother = GaussSeidel(),
    postsmoother = GaussSeidel(),
    kwargs...) where {T,V,bs,TA<:SparseMatrixCSC{T,V}}

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

        A = _extend_hierarchy!(levels, fine_fespace, coarse_fespace, A)

        coarse_x!(w, size(A, 1))
        coarse_b!(w, size(A, 1))
        residual!(w, size(A, 1))
    end
    return MultiLevel(levels, A, pcoarse_solver(A), presmoother, postsmoother, w)
end

function _extend_hierarchy!(levels, fine_fespace::FESpace, coarse_fespace::FESpace, A)
    P = build_prolongator(fine_fespace, coarse_fespace)
    R = P' # TODO: do we need other method to compute R?
    push!(levels, Level(A, P, R))
    A = R * A * P # Galerikn projection
    
    # rediscretization approach
    dropzeros!(A)
    return A
end
