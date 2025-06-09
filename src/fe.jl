## This file encompasses all the finite element utils required for P-Multigrid
struct FESpace{DH<:AbstractDofHandler,CV<:AbstractCellValues}
    dh::DH
    cv::CV
end

order(fe_space::FESpace) = fe_space.cv.fun_values.ip |> getorder
interpolation(fe_space::FESpace) = fe_space.cv.fun_values.ip
quadraturerule(fe_space::FESpace) = fe_space.cv.qr
Ferrite.ndofs(fe_space::FESpace) = ndofs(fe_space.dh)
Ferrite.getnbasefunctions(fe_space::FESpace) = getnbasefunctions(fe_space.cv)

function coarsen_order(fe_space::FESpace, p::Int)
    dh = fe_space.dh
    cv = fe_space.cv

    @assert 1 ≤ p < getorder(cv) "Invalid order $p for coarsening"

    # FIXME: more robust way to handle this?
    qr = fe_spqace |> quadraturerule
    ip = _new_coarse_ip(fe_space |> interpolation, p)
    coarse_cv = CellValues(qr, ip)
    coarse_dh = DofHandler(dh.grid)
    add!(coarse_dh, :x, ip) # TODO: do we need to use same name for the field?
    close!(coarse_dh)
    return FESpace(coarse_dh, coarse_cv)
end

function _new_coarse_ip(ip::Interpolation, p::Int)
    T = typeof(ip)
    BasisFunction = T.name.wrapper  #TODO: all basis functions have the same construction structure?
    RefShape = T.parameters[1]
    return BasisFunction{RefShape,p}()
end
