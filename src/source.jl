# Consider creating Source module in src.jl under src/ directory.  Such a module can be used
# by
#    include("src/src.jl")
#    importall (or using?) .Source
# See "importall .LinAlg" in julia-0.6/base/sysimg.jl.

# When we discretize a source, a source located in the interval that begins or ends with a
# ghost point is cumbersome.  In such a case, we have to assign current density at the ghost
# point, which is impossible because there is no entry in the current density array that
# is assigned to ghost points.
#
# Depending on the boundary condition, we have to handle such a case differently.  Suppose
# we are assiging a point source corresponding to the primal field parallel to the boundary.
# We assign a portion of the source current at the nonghost point adjacent to the ghost
# point.  The question is how to handle the portion of the source current to be assigned to
# the ghost point.
#
# - For Bloch, we have to assign a point source at the nonghost point corresponding to the
# ghost point.  Note that this should be additive to the existing current density at the
# nonghost point.  (The same situation as two point sources in neighboring grid intervals:
# at the shared grid point, the contributions from the two point sources are added.  The
# first and last grid intervals are adjacent through the Bloch boundary condition.)
#
# - For the symmetry boundary, putting current density component parallel to the boundary
# effectively puts another anti-symmetric current density component at the mirror symmetry
# position.  Usually these two current density components are sufficiently well-separated
# and their contributions do not overlap.  However, in the situation considered here, they are
# actually in neighboring cell intervals that share a point on the symmetry boundary.
#
# Now, when one point source is within a cell interval, the method to assign current density
# to the two end points of the cell interval is (a) to find the current density we would get
# when the source is exactly on one of the two points, and (b) calculate the portion of the
# total current to assign, depending on the location of the point in the cell interval: if
# the point is closer to one end point, assign a larger portion of the total current to that
# point and less to the other end point.
#
# As the two point sources with opposite directions approach to the symmetry boundary, the
# portions of the two point sources' currents assigned at the shared point on the boundary
# exactly cancel each other.  Therefore, we always assign 0 current density on the ghost
# point.  This makes us to simply ignore assiging current density to the point on the
# symmetry boundary.  The portion of the total current to assign to the grid point adjacent
# to the symmetry boundary is as usual, and this portion becomes zero as the point source
# approaches the symmetry boundary.  This makes sense, because if the point source is
# exactly on the symmetry boundary, being parallel to the boundary, then there must be no
# current in the system at all.
#
#
# Now, consider assignment of a point source corresponding to the dual field parallel to the
# boundary.  For Bloch, we proceed similarly.  For the symmetry boundary, now the imaginary
# source put behind the boundary is in the same direction as the original source in the
# domain.  Also, when the source is sufficiently close to the boundary, the actual and
# imaginary sources are in the same dual grid cell interval (whose end points are dual grid
# points), instead of the two adjacent cell intervals.  We don't have to assign a portion
# of the total current of the point source to the point behind the symmetry boundary
# ourselves: this is something that is automatically done by imposing symmetry.  However,
# if the point source existing behind the symmetry boundary assigns some portion of its
# total current to the point within the domain, that must be assigned ourselves.
#
# Because the two point sources are at mirror-symmetric locations, the reduction in the
# amount of the current assigned to the dual grid point within the domain by the point
# source moving away is exactly compensated by the increase in the amount of the current
# assigned to this point by the imaginary point source approaching to this point.
# Therefore, the net effect is the constant amount of current assigend to the dual grid
# point, no matter where the point source is located inside this cell interval containing
# the symmetry boundary.

# We need to support both Je and Jm to realize unidirectional sources.

export Source

# In a given Cartesian direction, if the source is located between two grid points,
# calculate the weights to distribute over the two points.  The output is in the form of
# ([ind₁, ind₂], [wt₁, wt₂]), where wtₖ is the weight assigned to the grid point at indₖ.
#
# If the source is placed exactly at a grid point, weight distribution is unnecessary and
# the output is in the form of ([ind], [wt]).  Note that the output type remains the same as
# before.
#
# This function works for both longitudinal and transversal discretization with respect to
# polarization.  Luckily, the above-described consideratino of the boundary conditions works
# for both polarization.  For example, suppose we are distributing electric current along
# the x-grid.  For Bloch, it is obvious that both polarizations are treated equally.  For
# the symmetry boundary, we can handle both polarization with the same code by determining
# whether the boundary condition is congruent with the polarization or not.  (See the
# definition of `withcongbc` in the function definition.)
#
# - For Jy and Jz, the passed grid locations are primal, so the congruent boundary is the
# symmetry boundary for primal fields.
# - For Jx, the passed grid locatinos are dual, so the congruent boundary is the symmetry
# boundary for dual fields.
#
# For both cases,
# - if the boundary is congruent, the source put beyond the boundary is antisymmetric to the
# source put before the boundary, so the two sources are superposed destructively at the
# boundary (note the two sources are in adjacent cells), and
# - if the boundary is incongruent, the source put beyond the boundary is symmetric to the
# source put before the boundary, so the two sources are superposed constructively at the
# boundary (note the two sources are in the same cell).
#
# In the implementation, assume that the source value to distribute is something that does
# not change with discretization (e.g., Id of a point source and K of a plane source).  The
# returned weights are the factors to multiply to this source value.
#
# When distributing a source between two points, the Bloch phase factor is not multiplied.
# The source we are trying to distribute by this function is a delta-function-like source
# in the axis of consideration (e.g., axis normal to a plane source).  It feels a bit
# strange to impose phase difference between what were originally a single point.  Also,
# I think it would be much easier to multiply the Bloch phase factors to the entire 3D array
# of source when all the sources are assigned and properly distributed with weight factors.
# That way, we don't have to implement multiplication of the Bloch phase factors in every
# source type.  However, this decision can change later if needed.
function distweights(c::Real,  # location of point (c means center)
                     domain::VecFloat,  # [lneg, lpos]: negative- and positive-end of the domain
                     l::VecFloat,  # location of grid planes without ghost points
                     ∆l::VecFloat,  # size of grid cell centered at l (not lg)
                     gt::GridType,  # type of the grid planes
                     isbloch::Bool)  # boundary condition
    length(domain)==2 || throw(ArgumentError("domain = $domain must be length-2."))
    domain[1]≤c≤domain[2] || throw(ArgumentError("c = $c must be within domain = $domain."))

    L = domain[2] - domain[1]
    N = length(l)

    # Below, what if N = 1?  I think N = 1 will work fine for Bloch, for which indn = indp
    # = 1 and therefore the contribution from indn and indp will be added in source array
    # construction, and for the symmetry boundary that is incongruent with the field type,
    # for which indn = indp = 1 again but only indn will be considered and the current
    # density to set will remain constant regardless of the source position... Wait, we have
    # contribution from both boundaries in this case, so maybe we should consider one real
    # point source and two imaginary point sources...
    #
    # Anyway, I think the main reason to use N = 1 is for 2D simulation, and in that case
    # we use Bloch in that direction.  So, I think it should be OK to allow N = 1 only for
    # Bloch.
    N > 1 || isbloch || throw(ArgumentError("length(l) = $N must be > 1 for symmetry boundary (= non-Bloch)."))

    # Test if c is in the range where the result is not affected by the boundary.
    # - For PRIM,
    #       - for Bloch, l[1] ≤ c < l[end], and
    #       - for symmetry, l[2] ≤ c < l[end].
    # - For DUAL,
    #       - for Bloch, l[1] ≤ c < l[end], and
    #       - for symmetry, l[1] ≤ c < l[end].
    # Therefore, only the case with the primal field and symmetry boundary is different.
    # I will call the a boundary condition congruent with the field type, because it is the
    # boundary that is effective for the field type of interest.
    cong_bc = gt==PRIM && !isbloch
    indn = cong_bc ? 2:1
    indp = N

    # Below, note that c < l[indp] is used rather than c ≤ l[indp].  If c = l[indp], then
    # ind₁ = findlast(x->x≤0, l.-c) is indp, and ind₂ = ind₁ + 1 is not well-defined when
    # indp = N.  The use of c < l[indp] avoids such a case.
    bc_noeff = l[indn] ≤ c < l[indp]  # boundary condition does not affect result
    if bc_noeff
        ind₁ = findlast(x->x≤0, l.-c)
        ind₂ = ind₁ + 1
        ∆c1 = c - l[ind₁]
        ∆lc = l[ind₂] - l[ind₁]
    else  # boundary affects result
        if c < l[indn]  # cannot happen for PRIM && Bloch
            ind₁ = indn
            ind₂ = N
            bnd = domain[1]
        else # c ≥ l[indp] (not >)
            ind₁ = indp
            ind₂ = 1
            bnd = domain[2]
        end

        ∆c1 = abs(l[ind₁] - c)
        if gt == PRIM
            ∆lc = abs(l[ind₁] - bnd)
        else  # gt == DUAL
            ∆lc = isbloch ? L-abs(l[ind₁]-l[ind₂]) : 2abs(l[ind₁]-bnd)
        end
    end

    r = ∆c1 / ∆lc

    # Below, scale wt₁ and wt₂ with 1-r and r, respectively,
    # - if c is in the region where the result is not affected by the boundary (bc_noeff), or
    # - even if c is in the region where the result is affected by the boundary,
    #       - if the boundary is Bloch, or
    #       - if the boundary is symmetry it is congruent with field type (cong_bc = true).
    wt₁ = 1.0 / ∆l[ind₁]
    (bc_noeff || isbloch || cong_bc) && (wt₁ *= 1.0-r)

    # Return only one weight factor
    # - if c is exactly at a grid point, or
    # - even if c is not exactly at a grid point, if it is in the region where the result is
    # affected by the boundary but the boundary is non-Bloch (= symmetry).
    if c==l[ind₁] || (!bc_noeff && !isbloch)
        # Return only one weight factor.
        ind = [ind₁]
        wt = [wt₁]
    else
        # Return two weight factors.  Make sure all the cases handled by this else block
        # have ind₂ defined already.
        wt₂ = 1.0 / ∆l[ind₂]
        (bc_noeff || isbloch) && (wt₂ *= r)

        ind = [ind₁, ind₂]
        wt = [wt₁, wt₂]
    end

    return (ind, wt)
end

abstract type Source end

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
    p::SVec3Float  # polarization direction
    c::SVec3Float  # location (c menas center)
    Id::Float  # dipole current I⋅d; default value = 1
end


mutable struct PlaneSrc <: Source
    n::Axis  # direction normal to plane
    ϕ::Float  # angle indicating polarization direction around n; ϕ = 0 is Cartesian direction cyclically next to n (e.g., for n = Ŷ, ϕ = 0 is z-direction)
    c::Float  # intercept (c means center)
    K::Float  # sheet current density (current per unit cross-sectional length in plane); default value = 1
end
