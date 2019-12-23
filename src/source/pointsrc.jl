export PointSrc

# When put a point source in space discretized by a finite-difference grid, the point source
# must specify current dipole I⋅d rather current I or current density J.
#
# - If it assigns current density J, it is assumed to flow uniformly through the area
# element ∆a where J is centered, so the total current is I = J ∆a, which decreases with
# increasing resolution.  If J is electric current density and if it decreases into a half,
# the induced magnetic field that circuits around I decreases into a half at the same
# distance from I.  Therefore, in order to keep the simulation result independent of grid
# resolution for sufficiently high resolution (i.e., in order to make the simulation result
# converges), current I instead of current density J should be assigned.
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
# In conclusion, a point source must specify the dipole current I⋅d, and from that it must
# recover J = I⋅d / (∆a ∆l) to assign to the field point.
#
# Wait, if we put a z-directional electric current component in a system that is periodic
# and one-cell-thick in the z-direction, shouldn't Jz be independent of the thickness ∆z of
# the system?  If we put a point source described above, then changing ∆z will change Jz to
# assign.  This is not what we want in general in simulating a 2D system.
#
# This means that a point source is not the correct source type to simulate a 2D system.
# This is obvious if we think about the actual physical systems represented by discretized
# systems.  The physical source distribution represented by the above-described situation is
# a periodic array of point sources separated by ∆z.  Here, by reducing ∆z while keeping the
# strength of each point source, we actually change the physical situation.  We have more
# charges per unit z oscillating along the z-direction; in other words, the charge density
# per unit z along the source line increases.  This change of charge density with varying ∆z
# should be captured in the discretized system, and the above-described point source
# specification in terms of the current dipole I⋅d does so: as ∆z decreases, Iz = Jz ∆x ∆y
# = I⋅d / ∆z increases indeed.
#
# On the other hand, the similar 2D situation we usually want to simulate is a 3D system
# with a z-directional line source with a constant current Iz.  Here, choosing ∆z is just to
# pick a certain ∆z-long section of the current, in which the current Iz does not change
# with the choice of ∆z.  Therefore, we must have a line source type that is specified with
# a current I instead of a current dipole I⋅d.  Using such a line source would not change
# Jz = I / (∆x ∆y) to assign on each edge element even if ∆z changes.

mutable struct PointSrc <: Source
    # Specify geometry first, current next (because the value of current is optional).
    c::SFloat{3}  # location (c menas center)
    p::SFloat{3}  # polarization direction
    Id::CFloat  # source strength (= dipole current I⋅d); default value = 1; complex in order to represent phase difference between point sources

    PointSrc(c::AbsVecReal, p::AbsVecReal, Id::Number=1.0) = new(c, normalize(p), Id)
end

PointSrc(c::AbsVecReal, p::Axis, Id::Number=1.0) = PointSrc(c, XYZ.==p, Id)

function add!(j3d::AbsArrNumber{4},  # 4D array of Je (electric current density) or Jm (magnetic current density)
              ft::FieldType,  # type of source (electric or magnetic)
              boundft::SVector{3,FieldType},  # boundary field type
              bounds::Tuple2{SFloat{3}},  # bounds[NEG][k] = boundary of domain at negative end in k-direction
              l::Tuple23{<:AbsVecFloat},  # l[PRIM][k] = primal vertex locations in k-direction
              ∆l::Tuple23{<:AbsVecFloat},  # ∆l[PRIM][k] = (∆l at primal vertices in w) == diff(l[DUAL][k] including ghost point)
              isbloch::SBool{3},  # Bloch boundary conditions
              src::PointSrc)  # point source to add
    gt_cmp₀ = ft2gt.(ft, boundft)  # grid type of source type (electric or magnetic)
    for nr = nXYZ  # set Jr
        Ird = src.Id * src.p[nr]  # r-component of I⋅d

        gt_cmp = broadcast((k,w,g)->(k==w ? alter(g) : g), nXYZ, nr, gt_cmp₀)  # grid type of Jr
        lcmp = t_ind(l, gt_cmp)
        ∆lcmp = t_ind(∆l, gt_cmp)

        x, ∆x = lcmp[nX], ∆lcmp[nX]
        indx, wtx = distweights(src.c[nX], gt_cmp[nX], t_ind(bounds,nX,nX), x, ∆x, isbloch[nX])

        y, ∆y = lcmp[nY], ∆lcmp[nY]
        indy, wty = distweights(src.c[nY], gt_cmp[nY], t_ind(bounds,nY,nY), y, ∆y, isbloch[nY])

        z, ∆z = lcmp[nZ], ∆lcmp[nZ]
        indz, wtz = distweights(src.c[nZ], gt_cmp[nZ], t_ind(bounds,nZ,nZ), z, ∆z, isbloch[nZ])

        wt_vxl = SArray{Tuple{2,1,1}}(wtx) .* SArray{Tuple{1,2,1}}(wty) .* SArray{Tuple{1,1,2}}(wtz)  # 2×2×2 SArray
        Jr = Ird .* wt_vxl

        for k = 1:2, j = 1:2, i = 1:2
            j3d[indx[i],indy[j],indz[k],nr] += Jr[i,j,k]
        end
    end

    return nothing
end
