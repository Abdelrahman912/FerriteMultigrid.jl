function pmultigrid(
    A::TA,
    fe_space::FESpace,
    ::Type{Val{bs}} = Val{1};
    presmoother = GaussSeidel(),
    postsmoother = GaussSeidel(),
    coarse_solver = AMGCoarseSolver, # TODO: remove Pinv and add AMGCoarseSolver
) where {T,V,bs,TA<:SparseMatrixCSC{T,V}}

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
    return MultiLevel(levels, A, coarse_solver(A), presmoother, postsmoother, w)
end

function _extend_hierarchy!(levels, fine_fespace::FESpace, coarse_fespace::FESpace, A)
    P = build_prolongator(fine_fespace, coarse_fespace)
    R = P' # TODO: do we need other method to compute R?
    push!(levels, Level(A, P, R))
    #A = R * A * P # Galerikn projection

    # rediscretization approach
    A = _rediscretize(coarse_fespace)
    #dropzeros!(A)
    return A
end


function _rediscretize(fe_space::FESpace)
    dh = fe_space.dh
    cellvalues = fe_space.cv
    grid = fe_space.dh.grid
    K = allocate_matrix(dh)

    ch = ConstraintHandler(dh)

    ∂Ω = union(
        getfacetset(grid, "left"),
        getfacetset(grid, "right")
    )

    dbc = Dirichlet(:x, ∂Ω, (x, t) -> 0.0)
    add!(ch, dbc)
    close!(ch)

    K, f = _assemble_global(cellvalues, K, dh)
    apply!(K, f, ch)
    return K
end

function _assemble_element!(Ke, fe, cellvalues)
    fill!(Ke, 0.0)
    fill!(fe, 0.0)
    n_basefuncs = getnbasefunctions(cellvalues)
    for q in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q)
        for i in 1:n_basefuncs
            ∇δu = shape_gradient(cellvalues, q, i)
            δu = shape_value(cellvalues, q, i)
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q, j)
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
            fe[i] += δu * dΩ  # RHS = 1
        end
    end
    return Ke, fe
end

function _assemble_global(cellvalues, K, dh)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        _assemble_element!(Ke, fe, cellvalues)
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end
