## 1D poisson equation with Dirichlet boundary conditions ##
function poisson(N, p, nqr)
    grid = generate_grid(Line, (N,))
    ip = Lagrange{RefLine, p}() #p2
    qr = QuadratureRule{RefLine}(nqr)  
    cellvalues = CellValues(qr, ip)

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    K = allocate_matrix(dh)

    ch = ConstraintHandler(dh)

    ∂Ω = union(
        getfacetset(grid, "left"),
        getfacetset(grid, "right")
    )

    dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0.0)
    add!(ch, dbc)
    close!(ch)

    K, f = _assemble_global(cellvalues, K, dh)
    apply!(K, f, ch)

    fe_space = FESpace(dh, cellvalues)
    return K, f, fe_space
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
