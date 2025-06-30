## defines how we obtain A operator for the coarse grid
abstract type AbstractCoarseningStrategy end
struct Galerkin <: AbstractCoarseningStrategy end
struct Rediscretization{TP <: AbstractPMultigrid} <: AbstractCoarseningStrategy
    problem::TP
end

## defines how we project from fine to coarse grid
abstract type AbstractProjectionStrategy end
struct DirectProjection <: AbstractProjectionStrategy end
struct StepProjection <: AbstractProjectionStrategy 
    step::Int
    function StepProjection(step::Int)
        step < 1 && error("Step must be greater than or equal to 1")
        return new(step)
    end
end

struct PMultigridConfiguration{TC<:AbstractCoarseningStrategy, TP<:AbstractProjectionStrategy}
    coarsening_strategy::TC
    projection_strategy::TP
end

pmultigrid_config() = PMultigridConfiguration(Galerkin(), DirectProjection())

function pmultigrid(
    A::TA,
    fe_space::FESpace,
    pgrid_config::PMultigridConfiguration ,
    pcoarse_solver , 
    ::Type{Val{bs}} = Val{1};
    presmoother = GaussSeidel(),
    postsmoother = GaussSeidel(),
    kwargs...) where {T,V,bs,TA<:SparseMatrixCSC{T,V}}

    levels = Vector{Level{TA,TA,Adjoint{T,TA}}}()
    w = MultiLevelWorkspace(Val{bs}, eltype(A))
    residual!(w, size(A, 1))
    
    p = fe_space |> order
    fespaces = Vector{FESpace}()
    push!(fespaces, fe_space)

    ps = pgrid_config.projection_strategy
    cs = pgrid_config.coarsening_strategy
    step  = _calculate_step(ps, p)

    while p > 1
        # reduce the polynomial order
        p = p - step > 1 ? p - step : 1

        fine_fespace = fespaces[end]
        coarse_fespace = coarsen_order(fine_fespace, p)
        push!(fespaces, coarse_fespace)

        A = _extend_hierarchy!(levels, fine_fespace, coarse_fespace, A, cs)

        coarse_x!(w, size(A, 1))
        coarse_b!(w, size(A, 1))
        residual!(w, size(A, 1))

       
    end
    return MultiLevel(levels, A, pcoarse_solver(A), presmoother, postsmoother, w)
end

function _extend_hierarchy!(levels, fine_fespace::FESpace, coarse_fespace::FESpace, A, ::Galerkin)
    P = build_prolongator(fine_fespace, coarse_fespace)
    R = P' # TODO: do we need other method to compute R?
    push!(levels, Level(A, P, R))
    A = R * A * P # Galerikn projection
    return A
end

function _extend_hierarchy!(levels, fine_fespace::FESpace, coarse_fespace::FESpace, A, cs::Rediscretization)
    P = build_prolongator(fine_fespace, coarse_fespace)
    R = P' # TODO: do we need other method to compute R?
    push!(levels, Level(A, P, R))

    problem = cs.problem
    A = assemble(problem, coarse_fespace)
    return A
end

function _calculate_step(ps::StepProjection, p::Int) 
    step = ps.step
    step ≥ p && error("Step must be less than the polynomial order $p")
    return step
end

_calculate_step(::DirectProjection, fine_p::Int) = fine_p - 1

