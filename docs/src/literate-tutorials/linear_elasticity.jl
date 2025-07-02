# # [Linear Elasticity](@id tutorial-linear-elasticity)
#
# ![](linear_elasticity.svg)
#
# *Figure 1*: Linear elastically deformed 1mm $\times$ 1mm Ferrite logo.
#
#md # !!! note
#md #     The full explanation for the underlying FEM theory in this example can be found in the [Linear Elasticity](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/linear_elasticity/) tutorial of the Ferrite.jl documentation.
#
#
# ## Implementation
# The following code is the same as the code in the [Linear Elasticity](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/linear_elasticity/) tutorial of the Ferrite.jl documentation, but with some comments removed for clarity.
# with two nuances:
# 1. Use of second order `Lagrange` shape functions for field approximation ⇒ `ip = Lagrange{RefTriangle,2}()^2`
# 2. Use 4 quadrature points to accomadate the second order shape functions ⇒ `qr = QuadratureRule{RefTriangle}(4)`
using Ferrite, FerriteGmsh, SparseArrays
using Downloads: download


function assemble_external_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)
    ## Create a temporary array for the facet's local contributions to the external force vector
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for facet in FacetIterator(dh, facetset)
        ## Update the facetvalues to the correct facet number
        reinit!(facetvalues, facet)
        ## Reset the temporary array for the next facet
        fill!(fe_ext, 0.0)
        ## Access the cell's coordinates
        cell_coordinates = getcoordinates(facet)
        for qp in 1:getnquadpoints(facetvalues)
            ## Calculate the global coordinate of the quadrature point.
            x = spatial_coordinate(facetvalues, qp, cell_coordinates)
            tₚ = prescribed_traction(x)
            ## Get the integration weight for the current quadrature point.
            dΓ = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                Nᵢ = shape_value(facetvalues, qp, i)
                fe_ext[i] += tₚ ⋅ Nᵢ * dΓ
            end
        end
        ## Add the local contributions to the correct indices in the global external force vector
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end

function assemble_cell!(ke, cellvalues, C)
    for q_point in 1:getnquadpoints(cellvalues)
        ## Get the integration weight for the quadrature point
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:getnbasefunctions(cellvalues)
            ## Gradient of the test function
            ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
            for j in 1:getnbasefunctions(cellvalues)
                ## Symmetric gradient of the trial function
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
            end
        end
    end
    return ke
end

function assemble_global!(K, dh, cellvalues, C)
    ## Allocate the element stiffness matrix
    n_basefuncs = getnbasefunctions(cellvalues)
    ke = zeros(n_basefuncs, n_basefuncs)
    ## Create an assembler
    assembler = start_assemble(K)
    ## Loop over all cells
    for cell in CellIterator(dh)
        ## Update the shape function gradients based on the cell coordinates
        reinit!(cellvalues, cell)
        ## Reset the element stiffness matrix
        fill!(ke, 0.0)
        ## Compute element contribution
        assemble_cell!(ke, cellvalues, C)
        ## Assemble ke into K
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end

function linear_elasticity_2d(C)
    logo_mesh = "logo.geo"
    asset_url = "https://raw.githubusercontent.com/Ferrite-FEM/Ferrite.jl/gh-pages/assets/"
    isfile(logo_mesh) || download(string(asset_url, logo_mesh), logo_mesh)

    grid = togrid(logo_mesh)
    addfacetset!(grid, "top", x -> x[2] ≈ 1.0) # facets for which x[2] ≈ 1.0 for all nodes
    addfacetset!(grid, "left", x -> abs(x[1]) < 1.0e-6)
    addfacetset!(grid, "bottom", x -> abs(x[2]) < 1.0e-6)

    dim = 2
    order = 2 # quadratic interpolation
    ip = Lagrange{RefTriangle,order}()^dim # vector valued interpolation

    qr = QuadratureRule{RefTriangle}(4) # 4 quadrature point
    qr_face = FacetQuadratureRule{RefTriangle}(1)

    cellvalues = CellValues(qr, ip)
    facetvalues = FacetValues(qr_face, ip)

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> 0.0, 2))
    add!(ch, Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> 0.0, 1))
    close!(ch)

    traction(x) = Vec(0.0, 20.0e3 * x[1])
    
    A = allocate_matrix(dh)
    assemble_global!(A, dh, cellvalues, C)

    b = zeros(ndofs(dh))
    assemble_external_forces!(b, dh, getfacetset(grid, "top"), facetvalues, traction)
    apply!(A, b, ch)

    return A, b, dh, cellvalues, ch
end


# ### Near Null Space
# Near Null Space (NNS) matrix for linear elasticity
function create_nns(dh)
    ##Ndof = ndofs(dh)
    grid = dh.grid
    Ndof = 2 * (grid.nodes |> length) # nns at p = 1 for AMG
    B = zeros(Float64, Ndof, 3)
    B[1:2:end, 1] .= 1 # x - translation
    B[2:2:end, 2] .= 1 # y - translation

    ## in-plane rotation (x,y) → (-y,x)
    coords = reduce(hcat, grid.nodes .|> (n -> n.x |> collect))' # convert nodes to 2d array
    y = coords[:, 2]
    x = coords[:, 1]
    B[1:2:end, 3] .= -y
    B[2:2:end, 3] .= x
    return B
end



using FerriteMultigrid
Emod = 200.0e3 # Young's modulus [MPa]
ν = 0.3        # Poisson's ratio [-]

Gmod = Emod / (2(1 + ν))  # Shear modulus
Kmod = Emod / (3(1 - 2ν)) # Bulk modulus

C = gradient(ϵ -> 2 * Gmod * dev(ϵ) + 3 * Kmod * vol(ϵ), zero(SymmetricTensor{2,2}))

A, b, dh, cellvalues, ch = linear_elasticity_2d(C);
B = create_nns(dh)
fe_space = FESpace(dh, cellvalues, ch)


## Galerkin Coarsening Strategy
config_gal = pmultigrid_config(coarse_strategy = Galerkin())
x_gal, res_gal = solve(K, f,fe_space, config_gal;B = B, log=true, rtol = 1e-10)

## Rediscretization Coarsening Strategy
config_red = pmultigrid_config(coarse_strategy = Rediscretization(LinearElasticityMultigrid(C)))
x_red, res_red = solve(K, f, fe_space, config_red; B = B, log=true, rtol = 1e-10)


using Test

@testset "Linear Elasticity Example" begin
    println("Final residual with Galerkin coarsening: ", res_gal[end])
    @test K * x_gal ≈ f

    println("Final residual with Rediscretization coarsening: ", res_red[end])
    @test K * x_red ≈ f
end



#md # ## [Plain program](@id linear-elasticity-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`linear_elasticity.jl`](linear_elasticity.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
