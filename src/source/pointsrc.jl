export PointSrc

# When put a point source in space discretized by a finite-difference grid, the point source
# must specify current dipole I d rather current I or current density J.
#
# - If it assigns current density J, it is assumed to flow uniformly through the area
# element ∆a where J is centered, so the total current is I = J ∆a, which decreases with
# increasing resolution.  If J is electric current and if it decreases into a half, the
# induced magnetic field that circuits around I decreases into a half at the same distance
# from I.  Therefore, in order to keep the simulation result independent of grid resolution
# for sufficiently high resolution (i.e., in order to make the simulation result converges),
# current I instead of current density J must be assigned.
#
# - Now, if it assigns current I, the current is assumed to flow uniformly in the edge
# element ∆l along the direction of the current.  (The fact that this assumption is made in
# finite difference is evident if we consider a loop of electric current formed along the
# perimeter of a pixel: the magnetic current dipole induced by this electric current loop is
# proportional to the line integral of the current, and on the finite-difference grid, this
# line integral is calculated by assuming each edge has a uniform current over its length.)
# For the same current I, if ∆l decreases, the current dipole I ∆l decreases.  When the
# current dipole decreases, the E-field induced by the dipole decreases.
#
# In conclusion, a point source must specify the dipole current I d, and from that it must
# recover J = I d / (∆a ∆l) to assign to the field point.

mutable struct PointSrc <: Source
    # Specify geometry first, current next (because the value of current is optional).
    c::SVec3Float  # location (c menas center)
    p::SVec3Float  # polarization direction
    Id::CFloat  # dipole current I⋅d; default value = 1

    PointSrc(c::AbsVecReal, p::AbsVecReal, Id::Number=1.0) = new(c, normalize(p), Id)
end

# function PointSrc(c::AbsVecReal, p::Axis, Id::Number=1.0)
#
#     ϕ = p==next1(n) ? 0.0 : π/2
#     return PlaneSrc(n, c, ϕ, K)
# end
#
# function add!(j3d::AbsArrNumber{4},  # 4D array of Je (electric current density) or Jm (magnetic current density)
#               gt::GridType,  # type of source (primal or dual)
#               bounds::Tuple2{SVec3Float},  # bounds[NEG][k] = boundary of domain at negative end in k-direction
#               l::Tuple23{<:AbsVecFloat},  # l[PRIM][k] = primal vertex locations in k-direction
#               ∆l::Tuple23{<:AbsVecFloat},  # ∆l[PRIM][k] = (∆l at primal vertices in w) == diff(l[DUAL][k] including ghost point)
#               isbloch::SVec3Bool,  # Bloch boundary conditions
#               src::PlaneSrc)  # sources; sources put later overwrite sources put earlier
#     # Set the coordinate axes.
#     nr = Int(src.n)
#     np, nq = next2(nr)
#
#     # Set the grid type and the staggered grid type
#     ngt = Int(gt)
#
#     # src.ϕ is the angle measured around r, src.ϕ = 0 corresponds to p.
#     Kp, Kq = cos(src.ϕ) * src.K, sin(src.ϕ) * src.K
#
#     lr, ∆lr = l[ngt][nr], ∆l[ngt][nr]
#     ind, wt = distweights(src.c, gt, t_ind(bounds,nr,nr), lr, ∆lr, isbloch[nr])
#     Jp, Jq = Kp*wt, Kq*wt  # Jp and Jq has one of two elements, like wt
#
#     jp3d, jq3d = @view(j3d[:,:,:,np]), @view(j3d[:,:,:,nq])
#     for k = 1:length(ind)
#         jp3d[Base.setindex(axes(jp3d), ind[k], nr)...] += Jp[k]
#         jq3d[Base.setindex(axes(jq3d), ind[k], nr)...] += Jq[k]
#     end
#
#     return nothing
# end
