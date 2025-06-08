using Ferrite
using AlgebraicMultigrid

using Ferrite, SparseArrays
using Test

# ---------------------------
# 1. Generate 1D grid (3 elements)
# ---------------------------
grid = generate_grid(Line, (3,))

# ---------------------------
# 2. Define quadratic Lagrange interpolation
# ---------------------------
ip = Lagrange{RefLine,2}()
qr = QuadratureRule{RefLine}(3)  # use enough quadrature points
cellvalues = CellValues(qr, ip)

# ---------------------------
# 3. Set up DofHandler
# ---------------------------
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# ---------------------------
# 4. Allocate matrix
# ---------------------------
K = allocate_matrix(dh)

# ---------------------------
# 5. Dirichlet BCs on both ends
# ---------------------------
ch = ConstraintHandler(dh)

∂Ω = union(getfacetset(grid, "left"), getfacetset(grid, "right"))

dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0.0)
add!(ch, dbc)
close!(ch)

# ---------------------------
# 6. Assemble element matrices
# ---------------------------
function assemble_element!(Ke, fe, cellvalues)
    fill!(Ke, 0.0)
    fill!(fe, 0.0)
    n_basefuncs = getnbasefunctions(cellvalues)
    for q = 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q)
        for i = 1:n_basefuncs
            ∇δu = shape_gradient(cellvalues, q, i)
            δu = shape_value(cellvalues, q, i)
            for j = 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q, j)
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
            fe[i] += δu * dΩ  # RHS = 1
        end
    end
    return Ke, fe
end

# ---------------------------
# 7. Assemble global matrix and RHS
# ---------------------------
function assemble_global(cellvalues, K, dh)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        assemble_element!(Ke, fe, cellvalues)
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

K, f = assemble_global(cellvalues, K, dh)

# ---------------------------
# 8. Apply constraints and solve
# ---------------------------
apply!(K, f, ch)
u = K \ f

# Prolongator
I = [1, 2, 2, 3, 4, 4, 5, 6, 6, 7];
J = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4];
#V = [1, 0.5, 0.5, 2, 0.5, 0.5, 2, 0.5, 0.5, 1];
V = [1, 0.5, 0.5, 2, 0.5, 0.5, 2, 0.5, 0.5, 1];
P = sparse(I, J, V)
R = P'

## test PolyMG
presmoother = GaussSeidel(iter = 4)
x = zeros(size(f))
presmoother(K, x, f)
r_h = f - K * x
r_2h = R * r_h
K_2h = R * K * P
e_2h = K_2h \ r_2h
e_h = P * e_2h
x += e_h
post_smoother = GaussSeidel(iter = 4)
post_smoother(K, x, f)
@test x ≈ u atol = 1e-1
