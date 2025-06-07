## M_ij = ∫ ϕ_i(x) ⋅ ϕ_j(x) dx
function mass_matrix!(Me::AbstractMatrix,cv::AbstractCellValues)
    fill!(Me, zero(eltype(Me)))
    n_basefuncs = getnbasefunctions(cv)
    for q in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q)
        for i in 1:n_basefuncs
            δu = shape_value(cv, q, i)
            for j in 1:n_basefuncs
                u = shape_value(cv, q, j)
                Me[i, j] += (δu * u) * dΩ   
            end
        end
    end
    return Me
end
