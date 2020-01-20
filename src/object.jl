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
export OpenInterval, ClosedInterval, KDTree, Object
export bounds, max∆l, matparam, paramind, objind, add!, periodize  #, surfpt_nearby, normal
# export lsf, bound_, L_, center_, dist2bound, bound_contains, ∆lmax, sphere, transform,
#     surfnormal, surfpoint  # functions
# import Base:size, getindex, contains, isless, union, intersect

mutable struct Object{K,Ke,Km,S<:Shape{K},Le,Lm}  # K: dimension of space, (Ke,Km): dimension of electric and magnetic material parameters
    shape::S
    mat::Material{Ke,Km,Le,Lm}
    ∆lmax::SVector{K,Float}
    oind::ObjInd  # object index; for comparison of sameness of objects quickly and more generally (e.g. when periodized; see add!)
    pind::Tuple2{ParamInd}  # {electric material index, magnetic material index} (see add!)
    Object{K,Ke,Km,S,Le,Lm}(shape, mat, ∆lmax) where {K,Ke,Km,S,Le,Lm} = new(shape, mat, ∆lmax)
end

Object(shape::S, mat::Material{Ke,Km}, ∆lmax::SVector{K}=@SVector(fill(Inf,K))) where {K,Ke,Km,S<:Shape{K}} =
    Object{K,Ke,Km,S,Ke*Ke,Km*Km}(shape, mat, ∆lmax)
Object(shape::Shape{K}, mat::Material, ∆lmax::AbsVec) where {K} = Object(shape, mat, SVector{K}(∆lmax))
Object(shape::Shape{K}, mat::Material, ∆lmax::Real) where {K} = Object(shape, mat, SVector(ntuple(k->∆lmax, Val(K))))

# Add a new convenience constructor
GeometryPrimitives.Box(b::Tuple2{AbsVec}, axes=Matrix{Float}(I,length(b[1]),length(b[1]))) = Box((b[1]+b[2])/2, abs.(b[2]-b[1]), axes)

GeometryPrimitives.bounds(o::Object) = bounds(o.shape)
Base.in(x::SVector{K}, o::Object{K}) where {K} = in(x, o.shape)
Base.in(x::AbsVec, o::Object{K}) where {K} = in(SVector{K}(x), o)

# Create the user interface such that maxwellsys.add(shape, material, ∆lmax) uses this function.
# setmat!(o::Object, m::Material) = (o.mat = m; o)
# setmax∆l!(o::Object{K}, ∆lmax::AbsVec) where {K} = (o.∆lmax = SVector{K}(∆lmax); o)
# setmax∆l!(o::Object{K}, ∆lmax::Number) where {K} = setmax∆l!(o, SVector(ntuple(k->∆lmax, Val(K))))

max∆l(o::Object) = o.∆lmax
matparam(o::Object, ft::FieldType) = matparam(o.mat, ft)
paramind(o::Object, ft::FieldType) = o.pind[Int(ft)]
objind(o::Object) = o.oind

# Define above functions for a vector of Object.  Define them with type parameters K, Ke, Km
# in order to make sure ovec is an array of Object with the same dimension.

# Even though Ke is not necessarily Km, the following function is still type-stable, because
# the output type is not a vector of SMatrix whose type depends on Ke or Km: the output is
# a regular array.
function matparam(ovec::AbsVec{Object{K,Ke,Km}}, ft::FieldType) where {K,Ke,Km}
    N = length(ovec)
    Kp = ft==EE ? Ke : Km
    paramvec = Array{CFloat,3}(undef, Kp, Kp, N)  # array representing vector of material parameter tensors

    for n = 1:N
        # It is debatable whether to index paramvec as [:,:,n] or [n,:,:].  paramvec is used
        # in assignment.jl.  The former is faster in iterating over the off-diagonal entries
        # inside assign_val!.  The latter is faster in iterating over a view of the (nw,nw)
        # diagonal entries as @view(paramvec[:,nw,nw]) because then the resulting vector
        # occupies a consecutive memory block.  I decide to use the former, because the
        # iteration over the off-diagonal entries is the innermost iteration occurring
        # inside assign_val!.  On the other hand, in setting the diagonal entries, once the
        # (nw,nw) entry of the parameter is taken, it is passed as a scalar to assign_val!
        # and iteration does not occur.
        paramvec[:,:,n] .= matparam(ovec[n], ft)
    end

    return paramvec
end

function paramind(ovec::AbsVec{Object{K,Ke,Km}}, ft::FieldType) where {K,Ke,Km}
    N = length(ovec)
    pindvec = Vector{ParamInd}(undef, N)

    for n = 1:N
        pindvec[n] = paramind(ovec[n], ft)
    end

    return pindvec
end

function objind(ovec::AbsVec{Object{K,Ke,Km}}) where {K,Ke,Km}
    N = length(ovec)
    oindvec = Vector{ObjInd}(undef, N)

    for n = 1:N
        oindvec[n] = objind(ovec[n])
    end

    return oindvec
end

# Consider using resize! on ovec.
function add!(ovec::AbsVec{Object{K,Ke,Km}}, paramset::Tuple2{AbsVec{SSComplex3}}, os::AbsVec{<:Object{K,Ke,Km}}) where {K,Ke,Km}
    for o = os
        add!(ovec, paramset, o)
    end
end

# Consider using resize! on ovec.
function add!(ovec::AbsVec{Object{K,Ke,Km}}, paramset::Tuple2{AbsVec{SSComplex3}}, os::Object{K,Ke,Km}...) where {K,Ke,Km}
    for o = os
        add!(ovec, paramset, o)
    end
end

# Add an object to the given vector `ovec` of objects and its electric and magnetic material
# parameters to the given sets `paramset` of material parameters.  In doing so, figure out
# the unique index for the object and the unique indices for the material parameters, and
# assign them to the object's fields `oind`  and `pind`.
#
# When I put an periodic array of an object, consider assigning the same object index to the
# periodized objects.  That way, I can treat two of objects over a periodic boundary as the
# same object.  (This is not implemented yet.)
#
# In other words, when assigning the same object index to distinct objects, we have to make
# sure that the two objects have the same surface normal directions at the equivalent points
# on the two objects, at least on the boundaries.  Therefore, for example if we put a cubic
# across a periodic boundary such that it is bisected by the boundary, we should not put one
# rectangular box that is tangential to one end of the domain and another rectangular box
# that is tengential to the other end of the domain.  Instead, we should create two cubics
# and place each of them centered on each boundary, such thay they stick out of the domain.
# (We can put two rectangular boxes that stick out of the domain, too, because they have to
# produce the same surface normals only at equivalent points on the boundaries.)
#
# Assigining the same object index is specifically to treat an object across a periodic
# boundary.  Therefore, we should not assign the same object index to the dintinct periodic
# objects in the domain (e.g., holes in a photonic crystal slab).
function add!(ovec::AbsVec{Object{K,Ke,Km}}, paramset::Tuple2{AbsVec{SSComplex3}}, o::Object{K,Ke,Km}) where {K,Ke,Km}
    # Assign the object index to o.
    # Currently, every object gets a new object index, but in the future different objects
    # may get the same object index; see the comments above about periodized objects.
    o.oind = isempty(ovec) ? 1 : objind(ovec[end])+1  # not just length(ovec)+1 to handle periodized objects in future
    push!(ovec, o)  # append o (for potential use of ovec with KDTree, must use pushfirst! to prepend)

    # Assign the material parameter indices to o.
    # If the material is aleady used in some object that was previously added, then use the
    # index already assigned to that material.
    p, pset = matparam(o,EE), paramset[nE]
    pind_e = findlast(x -> x==p, pset)  # number if p is already in pset
    pind_e==nothing && (push!(pset, p); pind_e = length(pset))  # update pind_e if p is not in pset

    p, pset = matparam(o,HH), paramset[nH]
    pind_h = findlast(x -> x==p, pset)  # number if p is already in psent
    pind_h==nothing && (push!(pset, p); pind_h = length(pset))  # update pind_h if p is not in pset

    o.pind = (pind_e, pind_h)

    return nothing
end

function GeometryPrimitives.periodize(o::Object{K}, A::AbsMat, ∆range::Shape{K}) where {K}
    shp_array = periodize(o.shape, A, ∆range)
    N = length(shp_array)
    obj_array = Vector{Object{K}}(N)
    for n = 1:N
        obj_array[n] = Object(shp_array[n], o.mat, o.∆lmax)
    end

    return obj_array
end


# The following includes operation on shapes.  These will need to be defined in GeometryPrimitives.
# Then, I will need to define the same functions for Object that delegate the operations to
# Shape's functions.
#
# See also test/object.jl


# # distance to the boundary.  If I want to get 0 for internal points, compare
# # dist(l,c) = (distance between l and center), and if it is less than L_(i) / 2,
# # return zero; otherwise return dist(l,c).
# @inline dist2bound(i::Interval1D, l::Real) = (a1 = abs(i.bound[nN]-l)) < (a2 = abs(i.bound[nP]-l)) ? a1 : a2

# # Functions common to all Shape's
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
