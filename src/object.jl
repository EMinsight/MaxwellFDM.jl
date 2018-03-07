# Below, the symbols of GeometryPrimitives such as KDTree, bounds, surfpt_nearby, normal do
# not have to be re-defined here to be exported.  When GeometryPrimitives is imported to
# MaxwellFDM by "using GeometryPrimitives", the symbols of GeometryPrimitives are already
# inside the namespace of MaxwellFDM, so we simply need to export them to use outside
# MaxwellFDM without qualification.
#
# One caveat, though.  In this file, bounds is defined for Interval, so if we simply export
# bounds, then only this is exported.  So, when defining bounds for Interval, we have to
# "extend" bounds for Shape by defining it as GeometryPrimitives.bounds(::Interval) = ....
# Then, exporting bounds exports this whole collection of bounds, both for Shape and Interval.
export OpenInterval, ClosedInterval, KDTree, Object, Object3
export bounds, max∆l, matparam, paramind, objind, add!  #, surfpt_nearby, normal
# export lsf, bound_, L_, center_, dist2bound, bound_contains, ∆lmax, sphere, transform,
#     surfnormal, surfpoint  # functions
# import Base:size, getindex, contains, isless, union, intersect

mutable struct Object{K,S<:Shape}  # use S<:Shape{K} after Julia issue #26321 is fixed
    shape::S
    mat::EncodedMaterial
    ∆lmax::SVector{K,Float}
    oind::Int  # object index; for comparison of sameness of objects quickly and more generally (e.g. when periodized; see add!)
    pind::Tuple2{Int}  # {primal material index, dual material index} (see add!)
    Object{K,S}(shape, mat, ∆lmax) where {K,S<:Shape{K}} = new(shape, mat, ∆lmax)
end

Object(shape::S, mat::EncodedMaterial, ∆lmax::SVector{K}=SVector(ntuple(k->Inf, Val{K}))) where {K,S<:Shape{K}} = Object{K,S}(shape, mat, ∆lmax)
Object(shape::Shape{K}, mat::EncodedMaterial, ∆lmax::AbsVec) where {K} = Object(shape, mat, SVector{K}(∆lmax))
Object(shape::Shape{K}, mat::EncodedMaterial, ∆lmax::Real) where {K} = Object(shape, mat, SVector(ntuple(k->∆lmax, Val{K})))

# Define Object with ∆lmax, Shape, and a material.
const Object3 = Object{3}

# Add a new convenience constructor
GeometryPrimitives.Box(b::Tuple2{AbsVec}, axes=eye(length(b[1]))) = Box((b[1]+b[2])/2, abs.(b[2]-b[1]), axes)

GeometryPrimitives.bounds(o::Object) = bounds(o.shape)
Base.in(x::SVector{K}, o::Object{K}) where {K} = in(x, o.shape)
Base.in(x::AbsVec, o::Object{K}) where {K} = in(SVector{K}(x), o)

# Create the user interface such that maxwellsys.add(shape, material, ∆lmax) uses this function.
# setmat!(o::Object, em::EncodedMaterial) = (o.mat = em; o)
# setmax∆l!(o::Object{K}, ∆lmax::AbsVec) where {K} = (o.∆lmax = SVector{K}(∆lmax); o)
# setmax∆l!(o::Object{K}, ∆lmax::Number) where {K} = setmax∆l!(o, SVector(ntuple(k->∆lmax, Val{K})))

max∆l(o::Object) = o.∆lmax
matparam(o::Object, gt::GridType) = matparam(o.mat, gt)
paramind(o::Object{K}, gt::GridType) where {K} = o.pind[Int(gt)]
objind(o::Object) = o.oind

function add!(ovec::AbsVec{<:Object{K}}, paramvec::Tuple2{AbsVec{SMat3Complex}}, os::Object{K}...) where {K}
    for o = os
        add!(ovec, paramvec, o)
    end
end

# When I put an periodic array of an object, consider assigning the same object index to the
# periodized objects.  That way, I can treat two of objects over a periodic boundary as the
# same object.
function add!(ovec::AbsVec{<:Object{K}}, paramvec::Tuple2{AbsVec{SMat3Complex}}, o::Object{K}) where {K}
    # Assign the object index to o.
    o.oind = isempty(ovec) ? 1 : objind(ovec[1])+1  # not just length(ovec)+1 to handle periodized objects
    unshift!(ovec, o)  # prepend o (for potential use of ovec with KDTree)

    # Assign the material parameter indices to o.
    p, pvec = matparam(o,PRIM), paramvec[nPR]
    pind_prim = findlast(x -> x==p, pvec)
    pind_prim==0 && (push!(pvec, p); pind_prim = length(pvec))

    p, pvec = matparam(o,DUAL), paramvec[nDL]
    pind_dual = findlast(x -> x==p, pvec)
    pind_dual==0 && (push!(pvec, p); pind_dual = length(pvec))

    o.pind = (pind_prim, pind_dual)

    return nothing
end


abstract type Interval end

struct OpenInterval <: Interval
    bounds::Tuple2{Float}
    ∆lmax::Float

    function OpenInterval(bounds::Tuple2{Real}, ∆lmax::Real=Inf)
        bounds[nP] ≥ bounds[nN] || throw(ArgumentError("bounds = $bounds must be ordered."))
        return new(bounds, ∆lmax)
    end
end
OpenInterval(o::Object, w::Int) = (b = bounds(o.shape); OpenInterval((b[nN][w], b[nP][w]), max∆l(o)[w]))  # w: direction

struct ClosedInterval <: Interval
    bounds::Tuple2{Float}
    ∆lmax::Float

    function ClosedInterval(bounds::Tuple2{Real}, ∆lmax::Real=Inf)
        bounds[nP] ≥ bounds[nN] || throw(ArgumentError("bounds = $bounds must be ordered."))
        return new(bounds, ∆lmax)
    end
end
ClosedInterval(o::Object, w::Int) = (b = bounds(o); ClosedInterval((b[nN][w], b[nP][w]), max∆l(o)[w]))  # w: direction

GeometryPrimitives.bounds(intv::Interval) = intv.bounds
max∆l(intv::Interval) = intv.∆lmax

Base.length(intv::Interval) = intv.bounds[nP] - intv.bounds[nN]
Base.in(l::Real, oi::OpenInterval) = oi.bounds[1] < l < oi.bounds[2]
Base.in(l::Real, ci::ClosedInterval) = ci.bounds[1] ≤ l ≤ ci.bounds[2]


# # distance to the boundary.  If I want to get 0 for internal points, compare
# # dist(l,c) = (distance between l and center), and if it is less than L_(i) / 2,
# # return zero; otherwise return dist(l,c).
# @inline dist2bound(i::Interval1D, l::Real) = (a1 = abs(i.bound[nN]-l)) < (a2 = abs(i.bound[nP]-l)) ? a1 : a2
#
# # Shapes
# abstract type Shape end
#
# immutable Box <: Shape
#     cbox::Interval3D  # circum-box
#     c::Tuple3{Float}  # center
#     s::Tuple3{Float}  # semiside
#
#     function Box(bound::Tuple32{Real}, ∆lmax::Tuple3{Real})
#         cbox = Interval3D(bound, ∆lmax)
#         c = center_(cbox)
#         s = L_(cbox) ./ 2
#
#         new(cbox, c, s)
#     end
# end
# Box(bound::Tuple32{Real}, ∆lmax::Real) = Box(bound, (∆lmax,∆lmax,∆lmax))
# Box(bound::Tuple32{Real}) = Box(bound, Inf)
#
# # Vectorization for a vector of l is achieved by the dot syntax: lsf.(b, l)
# lsf(b::Box) = (l::Tuple3{Real}) -> 1 - maximum(abs, (
#     (l[nX]-b.c[nX]) / b.s[nX],
#     (l[nY]-b.c[nY]) / b.s[nY],
#     (l[nZ]-b.c[nZ]) / b.s[nZ]
# ))
#
# immutable Ellipsoid <: Shape
#     cbox::Interval3D  # circum-box
#     c::Tuple3{Float}  # center
#     s::Tuple3{Float}  # semiaxes
#
#     function Ellipsoid(center::Tuple3{Real}, semiaxis::Tuple3{Real}, ∆lmax::Tuple3{Real})
#         any(semiaxis .≤ 0) && throw(ArgumentError("Elements of semiaxis = $semiaxis must be all positive."))
#
#         c = center
#         s = semiaxis
#         bound = ((c[nX]-s[nX],c[nX]+s[nX]), (c[nY]-s[nY],c[nY]+s[nY]), (c[nZ]-s[nZ],c[nZ]+s[nZ]))
#         cbox = Interval3D(bound, ∆lmax)
#
#         new(cbox, c, s)
#     end
# end
# Ellipsoid(center::Tuple3{Real}, semiaxis::Tuple3{Real}, ∆lmax::Real) = Ellipsoid(center, semiaxis, (∆lmax,∆lmax,∆lmax))
# Ellipsoid(center::Tuple3{Real}, semiaxis::Tuple3{Real}) = Ellipsoid(center, semiaxis, Inf)
#
# sphere(center::Tuple3{Real}, radius::Real, ∆lmax::Tuple3{Real}) = Ellipsoid(center, (radius,radius,radius), ∆lmax)
# sphere(center::Tuple3{Real}, radius::Real, ∆lmax::Real) = sphere(center, radius, (∆lmax,∆lmax,∆lmax))
# sphere(center::Tuple3{Real}, radius::Real) = sphere(center, radius, Inf)
#
# # Vectorization for a vector of l is achieved by the dot syntax: lsf.(b, l)
# lsf(e::Ellipsoid) = (l::Tuple3{Real}) -> 1 - hypot(
#     (l[nX]-e.c[nX]) / e.s[nX],
#     (l[nY]-e.c[nY]) / e.s[nY],
#     (l[nZ]-e.c[nZ]) / e.s[nZ]
# )
#
#
# # Functions common to all Shape's
#
# # size and getindex (to use the dot syntax for contains())
# size(::Shape) = ()
# getindex(s::Shape, ::CartesianIndex{0}) = s
#
# # contains
# # Sometimes, we only have level-set functions rather than Shape instances.
# @inline contains(lsf::Function, l::Tuple3{Real}, isinclusive::Bool=true) = isinclusive ? lsf(l) ≥ 0 : lsf(l) > 0
# @inline contains(s::Shape, l::Tuple3{Real}, isinclusive::Bool=true) = contains(lsf(s), l, isinclusive)
#
# # Other functions
# @inline bound_(s::Shape) = bound_(s.cbox)
# @inline bound_(s::Shape, w::Axis) = bound_(s.cbox, w)
# @inline bound_(s::Shape, w::Axis, sn::Sign) = bound_(s.cbox, w, sn)
# @inline bound_contains(s::Shape, l::Tuple3{Real}, isinclusive::Bool=true) = contains(s.cbox, l, isinclusive)
# @inline L_(s::Shape) = L_(s.cbox)  # scalar
# @inline center_(s::Shape) = center_(s.cbox)  # scalar
#
# # Level-set function manipulation
# @inline transform(lsf::Function, σ::Tuple3{Real}, ∆::Tuple3{Real}) = l -> lsf(σ .* (l .- ∆))  # σ[w] = 1 for shift, σ[w] = -1 for flip
#
# @inline union(lsf::Function...) = l -> maximum(map(f->f(l), lsf))
# @inline union(lsf1::Function, lsf2::Function) = lsf1==lsf2 ? lsf1 : l->max(lsf1(l), lsf2(l))
#
# @inline intersect(lsf::Function...) = l -> minimum(map(f->f(l), lsf))
#
# # Numerical gradient, and closest boundary point to a given point (needed for subpixel
# # smoothing of material parameters)
#
# # Calculates the direction normal to the curved surface crossing a voxel by numerical gradient
# # of the level-set function describing the surface.  The gradient is evaluated at a given point
# # inside the voxel.
# function surfnormal(lsf::Function, r::Tuple3{Float}, vxl::Tuple32{Float})
#     # println("lsf((0,0,0)) = $(lsf((0,0,0))),
#     #     lsf((-1,0,0)) = $(lsf((-1,0,0))), lsf((1,0,0)) = $(lsf((1,0,0))),
#     #     lsf((0,-1,0)) = $(lsf((0,-1,0))), lsf((0,1,0)) = $(lsf((0,1,0))),
#     #     lsf((0,0,-1)) = $(lsf((0,0,-1))), lsf((0,0,1)) = $(lsf((0,0,1))),
#     #     lsf((-1,-1,-1)) = $(lsf((-1,-1,-1))), lsf((1,-1,-1)) = $(lsf((1,-1,-1))),
#     #     lsf((-1,1,-1)) = $(lsf((-1,1,-1))), lsf((1,1,-1)) = $(lsf((1,1,-1))),
#     #     lsf((-1,-1,1)) = $(lsf((-1,-1,1))), lsf((1,-1,1)) = $(lsf((1,-1,1))),
#     #     lsf((-1,1,1)) = $(lsf((-1,1,1))), lsf((1,1,1)) = $(lsf((1,1,1))),
#     #     ")
#     lsf2(s) = 1-hypot(s[1], s[2], s[3])
#     # println("lsf($r) = $(lsf(r))")
#     # println("lsf2($r) = $(lsf2(r))")
#     # println("r = $r, vxl = $vxl")
#     α = 5e-4  # 1e-3 / 2
#
#     ∆x, ∆y, ∆z = α * (vxl[nX][nP] - vxl[nX][nN]), α * (vxl[nY][nP] - vxl[nY][nN]), α * (vxl[nZ][nP] - vxl[nZ][nN])
#     # ∆x, ∆y, ∆z = α .* (t_ind(vxl, nP, nP, nP) .- t_ind(vxl, nN, nN, nN))
#     rx, ry, rz = r
#
#     # println("(∆x, ∆y, ∆z) = $((∆x, ∆y, ∆z))")
#     # println("(rx, ry, rz) = $((rx, ry, rz))")
#     # # println("lsf(($rx, $ry, $rz)) = $(lsf((rx,ry,rz)))")
#     fxn::Float, fxp::Float = lsf((rx-∆x, ry, rz)), lsf((rx+∆x, ry, rz))
#     fyn::Float, fyp::Float = lsf((rx, ry-∆y, rz)), lsf((rx, ry+∆y, rz))
#     fzn::Float, fzp::Float = lsf((rx, ry, rz-∆z)), lsf((rx, ry, rz+∆z))
#     # println("((fxn, fxp), (fyn, fyp), (fzn, fzp)) = $(((fxn, fxp), (fyn, fyp), (fzn, fzp)))")
#
#     fxn2::Float, fxp2::Float = lsf2((rx-∆x, ry, rz)), lsf2((rx+∆x, ry, rz))
#     fyn2::Float, fyp2::Float = lsf2((rx, ry-∆y, rz)), lsf2((rx, ry+∆y, rz))
#     fzn2::Float, fzp2::Float = lsf2((rx, ry, rz-∆z)), lsf2((rx, ry, rz+∆z))
#     # println("((fxn2, fxp2), (fyn2, fyp2), (fzn2, fzp2)) = $(((fxn2, fxp2), (fyn2, fyp2), (fzn2, fzp2)))")
#
#     neg∂f∂x, neg∂f∂y, neg∂f∂z = (fxn - fxp) / ∆x, (fyn - fyp) / ∆y, (fzn - fzp) / ∆z
#
#     return neg∂f∂x, neg∂f∂y, neg∂f∂z  # outward normal is negative gradient
# end
#
# # Find the point on a curved surface by moving from a point along a direction.
# function surfpoint(lsf::Function, r::Tuple3{Float}, v::Tuple3{Float})
#     # println("r = $r, v = $v")
#     f(t) = lsf(r .+ t.*v)
#     tsol, isconverged = newtsol(0., f)
#
#     p = (r[nX] + tsol*v[nX], r[nY] + tsol*v[nY], r[nZ] + tsol*v[nZ])
#
#     return p, isconverged
#     # return r .+ tsol.*v
# end
