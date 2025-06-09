## P_ij = M^{-1}_{ij} ∫ ϕ_i(x) ⋅ ϕc_j(x) dx (ϕc is the coarse basis function)
function _element_prolongator!(
    Pe::AbstractMatrix,
    cv::AbstractCellValues,
    cv_coarse::AbstractCellValues,
)
    throw("Not implemented yet")
end

## M_{ij} = ∫ ϕ_i(x) ⋅ ϕ_j(x) dx
function _element_mass_matrix!(Me::AbstractMatrix, cv::AbstractCellValues)
    fill!(Me, zero(eltype(Me)))
    n_basefuncs = getnbasefunctions(cv)
    for q = 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q)
        for i = 1:n_basefuncs
            δu = shape_value(cv, q, i)
            for j = 1:n_basefuncs
                u = shape_value(cv, q, j)
                Me[i, j] += (δu * u) * dΩ
            end
        end
    end
    return Me
end

function _build_prolongator(fine_fespace::FESpace, coarse_fespace::FESpace)
    throw("Not implemented yet")
end

function build_prolongator(fine_fespace::FESpace, coarse_fespace::FESpace)
    error("Not implemented yet")
    # return P
end
