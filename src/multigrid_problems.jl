# interface for multigrid problems
# there will be ready to use implementations for some common problems,
# however for new problems, this interface should be implemented
abstract type AbstractPMultigrid end

function assemble(problem::AbstractPMultigrid, fe_space::FESpace)
    # this is an interface that should be defined for each specific problem type
    # This function should be implemented in the specific problem type
    error("assemble not implemented for $(typeof(problem))")
end

#######################
## Diffusion problem ##
#######################

abstract type AbstractCoefficient end
struct ConstantCoefficient{Tv <: Real} <: AbstractCoefficient
    K::Tv
end

function *(c::ConstantCoefficient, x::Real)
    return c.K * x
end

function *(x::Real, c::ConstantCoefficient)
    return x * c.K
end

struct DiffusionMultigrid{C} <: AbstractPMultigrid
    coeff::C
end

function DiffusionMultigrid(coeff::Real)
    return DiffusionMultigrid(ConstantCoefficient(coeff))
end


function assemble(problem::DiffusionMultigrid, fe_space::FESpace)
    dh = fe_space.dh
    cv = fe_space.cv
    ch = fe_space.ch

    K = _assemble_global_K(cv, dh, problem)
    apply!(K, ch)
    return K
end


function _assemble_element_K!(Ke, cellvalues,problem::DiffusionMultigrid)
    fill!(Ke, 0.0)
    n_basefuncs = getnbasefunctions(cellvalues)
    for q in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q)
        for i in 1:n_basefuncs
            ∇δu = shape_gradient(cellvalues, q, i)
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q, j)
                Ke[i, j] += problem.coeff * (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
end

function _assemble_global_K(cellvalues, dh,problem::DiffusionMultigrid)
    K = allocate_matrix(dh)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        _assemble_element_K!(Ke, cellvalues, problem)
        assemble!(assembler, celldofs(cell), Ke)
    end
    return K
end

