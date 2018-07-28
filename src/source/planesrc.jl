export PlaneSrc

mutable struct PlaneSrc <: Source
    # Specify geometry first, current next (because the value of current is optional).
    n::Axis  # axis normal to plane
    c::Float  # intercept (c means center)
    ϕ::Float  # angle indicating polarization direction around n; ϕ = 0 is Cartesian direction cyclically next to n (e.g., for n = Ŷ, ϕ = 0 is z-direction)
    K::CFloat # sheet current density (current per unit in-plane length normal to flow direction); default value = 1

    PlaneSrc(n::Axis, c::Real, ϕ::Real, K::Number=1.0) = new(n, c, ϕ, K)
end

function PlaneSrc(n::Axis, c::Real, p::Axis, K::Number=1.0)
    p≠n || throw(ArgumentError("p = $p must be orthogonal to n = $n."))
    ϕ = p==next1(n) ? 0.0 : π/2
    return PlaneSrc(n, c, ϕ, K)
end

function add!(j3d::AbsArrNumber{4},  # 4D array of Je (electric current density) or Jm (magnetic current density)
              gt::GridType,  # type of source (primal or dual)
              bounds::Tuple2{SVec3Float},  # bounds[NEG][k] = boundary of domain at negative end in k-direction
              l::Tuple23{<:AbsVecFloat},  # l[PRIM][k] = primal vertex locations in k-direction
              ∆l::Tuple23{<:AbsVecFloat},  # ∆l[PRIM][k] = (∆l at primal vertices in w) == diff(l[DUAL][k] including ghost point)
              isbloch::SVec3Bool,  # Bloch boundary conditions
              src::PlaneSrc)  # sources; sources put later overwrite sources put earlier
    # Set the coordinate axes.
    nr = Int(src.n)
    np, nq = next2(nr)

    # Set the grid type and the staggered grid type
    ngt = Int(gt)

    # src.ϕ is the angle measured around r, src.ϕ = 0 corresponds to p.
    Kp, Kq = cos(src.ϕ) * src.K, sin(src.ϕ) * src.K

    lr, ∆lr = l[ngt][nr], ∆l[ngt][nr]
    ind, wt = distweights(src.c, gt, t_ind(bounds,nr,nr), lr, ∆lr, isbloch[nr])
    Jp, Jq = Kp*wt, Kq*wt  # Jp and Jq has one of two elements, like wt

    jp3d, jq3d = @view(j3d[:,:,:,np]), @view(j3d[:,:,:,nq])
    for k = 1:length(ind)
        jp3d[Base.setindex(indices(jp3d), ind[k], nr)...] += Jp[k]
        jq3d[Base.setindex(indices(jq3d), ind[k], nr)...] += Jq[k]
    end

    return nothing
end
