# Consider creating Source module in source.jl under source/ directory.  Such a module can
# be used by
#    include("source/source.jl")
#    importall (or using?) .Source
# See "importall .LinAlg" in julia-0.6/base/sysimg.jl.

# When we discretize a source, a source located in the interval that begins or ends with a
# ghost point is cumbersome.  In such a case, we have to assign current density at the ghost
# point, which is impossible because there is no entry in the current density array that
# is assigned to ghost points.
#
# Depending on the boundary condition, we have to handle such a case differently.  Suppose
# we are assiging a point source of the primal field type in the direction parallel to (but
# slightly away from) the boundary.  We assign a portion of the source current at the
# nonghost point adjacent to the ghost point.  The question is how to handle the portion of
# the source current to be assigned to the ghost point.
#
# - For Bloch, we have to assign a point source at the nonghost point corresponding to the
# ghost point.  Note that this should be additive to any existing current density at the
# nonghost point.  (The same situation as two point sources in neighboring grid intervals:
# at the shared grid point, the contributions from the two point sources are added.  The
# first and last grid intervals are adjacent through the Bloch boundary condition.)
#
# - For the symmetry boundary, putting current density component parallel to the boundary
# effectively puts another anti-symmetric current density component at the mirror symmetry
# position.  Usually these two current density components are sufficiently well-separated
# and their contributions do not overlap.  However, in the situation considered here, they
# are actually in neighboring cell intervals that share a point on the symmetry boundary.
#
# If a point source is in a cell interval between two grid points A and B, the method to
# assign current densities to A and B is the "sliding" method.  In this method, imagine the
# point source location changes gradually from A to B.  When it is exactly at A, current
# density is assigned only to A, and the magnitude J_A of the current density is determined
# by dividing the current of the point source by the grid cell area occupied by A.  We can
# similarly calculate the current density J_B to assign to B in the case where the point
# source is exactly at B.  (Note that J_A ≠ J_B in general even though the current of the
# point source is the same, because the area occupied by A and B can be different if the
# grid is nonuniform.)  Now, if the point source location changes gradually from A to B, the
# current density to assign to A should decrease linearly from J_A to 0, and the current
# density to assign to B should increase linearly from 0 to J_B.  Therefore, if the point
# source is located between A and B such that the ratio of the distances to A and B is
# x:(1-x), then the current densities to assign to A and B should be J_A * (1-x) and J_B * x,
# respectively.  Note that these becomes J_A and 0 for x = 0 (when the point source is at A)
# and 0 and J_B for x = 1 (when the point source is at B).
#
# The situation of the symmetry boundary condition discussed earlier, where the two point
# sources with opposite directions are placed at the mirror symmetry locations around the
# symmetry boundary, can be thought of as assigning the current densities for the two
# individual point sources using the above described recipe.  Note that as the two point
# sources approach to the symmetry boundary, the portions of the two point sources' currents
# assigned at the shared point on the boundary exactly cancel each other.  Therefore, we
# always assign zero current density on the ghost point, which is at the symmetry boundary.
# This makes us to simply ignore assiging current density to the point on the symmetry
# boundary.  The portion of the total current to assign to the grid point adjacent
# to the symmetry boundary is as usual, and this portion becomes zero as the point source
# approaches the symmetry boundary.  This makes sense, because if the point source is
# exactly on the symmetry boundary, being parallel to the boundary, then there must be no
# current in the system at all.
#
#
#
# Now, consider assignment of a point source corresponding to the dual field parallel to the
# boundary.  For Bloch, we proceed similarly.  For the symmetry boundary, now the image
# source put behind the boundary is in the same direction as the original source in the
# domain.  Also, when the source is sufficiently close to the boundary, the actual and
# image sources are in a single dual grid cell interval, whose  the middle point is the
# boundary point and the end points are dual grid points, which are referred to as A and B
# as before for convenience below, instead of the two adjacent cell intervals.  Suppose A is
# within the domain and B is outside the domain.  We don't have to assign a portion of the
# total current of the real point source to B ourselves, which is behind the symmetry
# boundary: this is something automatically done by imposing symmetry.  However, if the
# image point source existing behind the symmetry boundary assigns some portion of its total
# current to A, which is within the domain, that must be assigned ourselves.
#
# Because the two point sources are at mirror-symmetric locations, the reduction in the
# amount of the current assigned to A by the point source moving away from A towards the
# symmetry boundary is exactly compensated by the increase in the amount of the current
# assigned to A by the image point source to the symmetry boundary from behind the boundary.
# Therefore, the net effect is the constant amount of current assigend to the dual grid
# point, no matter where the point source is located inside this cell interval containing
# the symmetry boundary.


# Note that we need to support both Je and Jm to realize unidirectional sources.

export Source
export add_src!

abstract type Source{K,Kf} end  # K: shape dimension; Kf: field dimension

add_src!(jKd::AbsArrNumber, gt₀::AbsVec{GridType}, bounds::Tuple2{AbsVecReal},
         l::Tuple2{NTuple{K,AbsVecReal}}, ∆l::Tuple2{NTuple{K,AbsVecReal}},
         isbloch::AbsVecBool, srcs::Source...) where {K} =
    add_src!(jKd, SVec{K}(gt₀), (float.(l[nPR]),float.(l[nDL])), (float.(∆l[nPR]), float.(∆l[nDL])), SVec{K}(isbloch), srcs...)

# About the determination of gt_cmp
# In the concrete implementation of add_src!() for each concrete Source type, gt_cmp, which
# is the grid type of the J field component to set up (typically indicated by the
# nw-component) is determined by
#
#     gt_cmp = src.isfield˔shp ? gt₀ : gt_w(nw, gt₀)
#
# gt_cmp is used to obtain the locations of the field by
#
#     lcmp = t_ind(l, gt_cmp)
#
# So, why do we choose gt₀ when the field subspace is orthogonal to the shape subspace?
#
# gt₀ represents the grid type of the corners of Yee's voxel whose edges are composed of the
# field lines.  The w-component of the field is located a half grid point away from these
# voxel corners along the w-direction.  This picture describes how to determine gt_cmp when
# the field subspace is equal to the shape subspace.
#
# However, when the field subspace is orthogonal to the shape subspace, the half grid point
# shift to apply to the w-component occurs in the direction normal to the shape subspace.
# Along that direction, the grid lines are not defined, because l is defined in the shape
# subspace.  Therefore, the half grid point shift does not have any effect in determining
# the locations of the w-component, and we can use gt₀ as the grid type in determining the
# field locations.
#
# For example, In the TE equation, the field subspace for the H-field is the z-direction and
# the shape subspace is the xy-plane.  When the boundary field types are the E-field type
# for both the x- and y-directions, gt₀ for the H-field (Hz) is [DUAL,DUAL].  We can use
# [DUAL,DUAL] in determining the xy-locations of the H-field, as the half grid point shift
# along the z-direction for the H-field does not affect the locations of the H-field.
function add_src!(jKd::AbsArrNumber{K₊₁},  # (K+1)-dimensional array of Je (electric) or Jm (magnetic); first K dimensions specify location; last 1 dimension specify field component
                  gt₀::SVec{K,GridType},  # grid type of voxel corners; generated by ft2gt.(ft, boundft)
                  bounds::Tuple2{SFloat{K}},  # bounds[NEG][k] = boundary of domain at negative end in k-direction
                  l::Tuple2{NTuple{K,VecFloat}},  # l[PRIM][k] = primal vertex locations in k-direction
                  ∆l::Tuple2{NTuple{K,VecFloat}},  # ∆l[PRIM][k] = (∆l at primal vertices in w) == diff(l[DUAL][k] including ghost point)
                  isbloch::SBool{K},  # Bloch boundary conditions
                  src::Source{K,Kf}...  # sources
                  ) where {K,K₊₁,Kf}
    for src = srcs
        add_src!(j3d, ft, boundft, bounds, l, ∆l, isbloch, src)
    end
end

# I could have defined distweights with more abstractly typed arguments, such as c::Real,
# bounds::AbsVecReal, l::AbsVecReal, ∆l::AbsVecReal, and have let the type inference system
# create distweights for each distinct set of types.  However, then I cannot control the
# behavior of the function perfectly, because the two numbers used in division could be
# integer or float.  Therefore, I implement the kernel for a specific type later, and here I
# define a wrapper that converts the arguments into those specific type.
#
# Another advantage of this approach is that even though several different versions of the
# wrapper will be generated, the kernel will be generated only once for the specific type.
# This makes the amount of compiled code minimal.
#
# Note that the situation here is a bit different from create_curl.  There, I wanted to
# create a curl operator of differnt types, such as Float, CFloat, and even Int.  Therefore,
# I had to allow generation of different kernel code for different argument types.
distweights(c::Real, gt::GridType, bounds::AbsVecReal, l::AbsVecReal, ∆l::AbsVecReal, isbloch::Bool) =
    distweights(float(c), gt, SFloat{2}(bounds), float(l), float(∆l), isbloch)

# In a given Cartesian direction, if the source is located between two grid points,
# calculate the weights to distribute over the two points.  The output is in the form of
# ([ind₁, ind₂], [wt₁, wt₂]), where wtₖ is the weight assigned to the grid point at indₖ.
# Note that the physical dimension of the weight is 1/(length).
#
# If the source is placed exactly at one grid point, weight distribution is unnecessary and
# the output satisfies ind₂ = ind₁ and wt₂ = 0.  Note that the output type remains the same
# as before, which is a useful feature for achieving type stability of the function calling
# this function.
#
# This function works for both longitudinal and transversal discretization with respect to
# polarization.  Luckily, the consideration of the boundary conditions discussed at the
# beginning of this file works for both cases.  For example, suppose we are distributing
# along the x-grid an electric current source.  For Bloch, it is obvious that both
# polarizations are treated equally.  For the symmetry boundary condition, we can handle
# both polarization with the same code by determining whether the boundary condition is
# "zeroing" the source on the boundary or not.  Specifically, the zeroing boundary condition
# is defined as follows (see the definition of `zeroing_bc` in the function definition as
# well):
# - if the boundary condition is zeroing the source, the source put beyond the boundary is
# antisymmetric to the source put within the boundary, so the two sources are superposed
# destructively on the boundary (note the two sources are in adjacent cells), and
# - if the boundary condition is not zeroing the source, the source put beyond the boundary
# is symmetric to the source put within the boundary, so the two sources are superposed
# constructively on the boundary (note the two sources are in the same cell).
#
# Therefore,
# - For Jy and Jz (parallel to the x-boundary), the zeroing boundary condition is the
# symmetric E-field boundary condition.  Because the x-location of Jy and Jz is the E-field
# grid location, we know that a symmetric boundary is the E-field boundary if the x-location
# of Jy and Jz are primal.
# - For Jx (normal to the x-boundary), the zeroing boundary condition is the symmetric
# H-field boundary condition.  Because the x-location of Jx is the H-field grid location, we
# know again that a symmetric boundary is the H-field boundary if the x-location of Jx is
# primal.
#
# Note that regardless of the polarization of the source, the boundary condition is zeroing
# if the location type is primal.
#
# In the implementation, assume that the source "value" to distribute is something that does
# not change with discretization (e.g., Id of a point source and K of a plane source).  The
# returned weights are the factors to multiply to this source value.
#
# When distributing a source between two points on the opposite sides of the Bloch boundary,
# the Bloch phase factor is not multiplied. The source we are trying to distribute by this
# function is a delta-function-like source in the axis of consideration (e.g., axis normal
# to a plane source).  It feels a bit strange to impose phase difference between what were
# originally a single point.  Also, I think it would be much easier to multiply the Bloch
# phase factors to the entire 3D array of source when all the sources are assigned and
# properly distributed with weight factors.  That way, we don't have to implement
# multiplication of the Bloch phase factors in every source type.  However, this decision
# can change later if needed.
function distweights(c::Float,  # location of point (c means center)
                     gt::GridType,  # type of the grid planes
                     bounds::SFloat{2},  # [lneg, lpos]: negative- and positive-end boundaries of domain
                     l::AbsVecFloat,  # location of grid planes without ghost points
                     ∆l::AbsVecFloat,  # size of grid cell centered at l (not lg)
                     isbloch::Bool)  # boundary condition
    bounds[1]≤c≤bounds[2] || throw(ArgumentError("c = $c must be within bounds = $bounds."))

    L = bounds[2] - bounds[1]
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
    # - For l[1:end] being the primal grid locations except for the ghost point,
    #       - for Bloch, l[1] ≤ c < l[end], and
    #       - for symmetry, l[2] ≤ c < l[end].  (Note the minimum boundary beng l[2], not l[1].)
    # - For l[1:end] being the dual grid locations except for the ghost point,
    #       - for Bloch, l[1] ≤ c < l[end], and
    #       - for symmetry, l[1] ≤ c < l[end].
    # Therefore, only the case with the primal field type and symmetry boundary condition is
    # different from other cases.  I will call this a case with a boundary condition
    # "congruent" with the field type, because it is the boundary that is effective for the
    # field type of interest.
    zeroing_bc = gt==PRIM && !isbloch
    indn = zeroing_bc ? 2 : 1
    indp = N

    # Below, note that c < l[indp] is used rather than c ≤ l[indp].  If c = l[indp], then
    # ind₁ = findlast(x->x≤0, l.-c) is indp, and ind₂ = ind₁ + 1 is not well-defined when
    # indp = N.  The use of c < l[indp] avoids such a case.
    bc_!eff = l[indn] ≤ c < l[indp]  # boundary condition does not affect result
    if bc_!eff
        ind₁ = findlast(x->x≤0, l.-c)  # because l is sorted and c ≥ l[1], l[1]-c ≤ 0 and thus ind₁ ≠ nothing
        ind₂ = ind₁ + 1
        ∆c1 = c - l[ind₁]
        ∆lc = l[ind₂] - l[ind₁]
    else  # boundary affects result
        if c < l[indn]  # cannot happen for PRIM && Bloch
            ind₁ = indn
            ind₂ = N
            bnd = bounds[1]
        else # c ≥ l[indp] (not >)
            ind₁ = indp
            ind₂ = 1
            bnd = bounds[2]
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
    # - if c is in the region where the result is not affected by the boundary (bc_!eff), or
    # - even if c is in the region where the result is affected by the boundary,
    #       - if the boundary is Bloch, or
    #       - if the boundary is symmetry it is congruent with field type (zeroing_bc = true).
    wt₁ = 1.0 / ∆l[ind₁]
    (bc_!eff || isbloch || zeroing_bc) && (wt₁ *= 1.0-r)

    # Return only one weight factor
    # - if c is exactly at a grid point, or
    # - even if c is not exactly at a grid point, if it is in the region where the result is
    # affected by the boundary but the boundary is non-Bloch (= symmetry).
    if c==l[ind₁] || (!bc_!eff && !isbloch)
        # We want to return only one weight factor in this case, but in order to keep the
        # output SVec type the same (SVec{2} rather than SVec{1}), introduce a
        # dummy second weight factor.
        ind₂ = ind₁
        wt₂ = 0.0
    else
        # Return two weight factors.  Make sure all the cases handled by this else block
        # have ind₂ defined already.
        wt₂ = 1.0 / ∆l[ind₂]
        (bc_!eff || isbloch) && (wt₂ *= r)
    end

    ind = SVec(ind₁, ind₂)
    wt = SVec(wt₁, wt₂)

    return (ind, wt)
end

# funciton calc_blochphase()
# function j3d2j()

# When defining concrete source types below, try to use the following order of arguments in
# the constructor:
# - geometry (e.g., direction normal to a plane)
# - location (e.g., location for a point, intercept for a plane)
# - polarization
# - strength (e.g., current dipole, current, sheet current density)

# I think I can define a "uniform source" that populates a uniform current inside a given
# geometrical object specified with a GeometryPrimitives type.  This can cover various
# primitive source types, such as a point source, line source, and a plane source.  The
# source assignment function will accept a vector of sources, assign these sources similarly
# as assign_param! in assignment.jl (i.e., iterate over a vector of shapes and assign points
# included in each shape).  A main difference is that I have to find voxels that overlap with
# the shape, instead of points in the shape.  Once such voxels are found, I need to assign
# source current density to the edges of the voxels.  (I have to think about subvoxel
# smoothing in a voxel that is partially included in the shape.)
#
# Even more generally, I think I can create a single source that accepts a geometry and
# source strength function.  The uniform sources are a special case with a constant
# source strength.

include("pointsrc.jl")
include("planesrc.jl")
