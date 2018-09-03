# Below, obj3d and oind3d are different: obj3d stores references to Object3 instances,
# whereas oind3d stores object indices.  The reason for this is to handle multiple objects
# forming essentially a single object by being connected through a periodic boundary.  In
# that case, we assign the same object index to those distinct objects.
#
# I could enhance the assignment performance by constructing oid3d (not oind3d), which stores
# unique object IDs rather than the reference itself.  This array is different from oind3d
# in that it distinguishes objects repeated by periodic boundary condition.  Then, I can
# retrieve objects from a map from this IDs to objects.  I could use oid3d in smoothing as
# well, by constructing an 8-vector oid3d_vxl, like pind3d_vxl and oind3d_vxl.
#
# The rationale for this trick is that assigning an Int to VecInt is faster than assigning a
# concrete Object (like Box) to Vector{Object3}.  In my test,
#   setindex!(::VecInt, ::Int, <index>)
# was about twice as fast as
#   setindex!(::Vector{Object3}, ::Box, <index>), though both took only a few nanoseconds.
#
# Currently the performance of assignment and smoothing seems satisfactory, so I will not
# pursue this extra optimization.  If I need to assign many objects, the situation may
# change and I may need to implement this optimization.

export create_param3d, create_n3d, assign_param!, assign_val_shape!
# using BenchmarkTools

# About the order of indices of param3d:
#
# param3d has 5 indices: the first three are positional indices (i,j,k) and the last two are
# material parameter tensor component indices (v,w).  This mean param3d[i,j,k,:,:] is the
# 3×3 material tensor.

# This is different from the conventional indexing scheme used for material parametetrs.
# For example, we usually write the xy-component of ε at a location (i,j,k) as ε_xy[i,j,k],
# where the tesnor component subscripts v = x and w = y comes earlier than the positional
# indices (i,j,k).  Similarly, the x-component of the E-field is usually written E_x[i,j,k],
# not E[i,j,k]_x.
#
# However, when assigning the material parameters in the code in this file, we ofter assign
# the same value to a block of param3d.  This block is chosen differently for different
# material tensor components, because different tensor components are evaluated at different
# locations in Yee's grid.  Therefore, in this assignment, we first fix the material tensor
# components v and w, determine the range of (i,j,k) to which the same material parameter
# value to assign, and perform the assignment.
#
# Now, if we index param3d as param3d[v,w,i,j,k], v and w are the fastest-varing indices.
# Therefore, fixing them and assigning for a contiguous range of (i,j,k) does not actually
# assign in a contiguous memory block.  On the other hand, if we index param3d as
# param3d[i,j,k,v,w], for fixed v and w a contiguous range of (i,j,k) actually assigns in
# a more contiguous memory block.  (The entries with adjacent i for the same j,k,v,w are
# actually adjacent in the memory space.)  Therefore, the latter indexing scheme results in
# faster performance.  In my experiment of assigning 8 objects in 3D, the latter turns out
# to be about 30% faster.
#
# For this reason, I use the less conventional indexing scheme of param3d[i,j,k,v,w].

# Creates param3d of size N+1.  It is created with size N+1 to match the sizes of other
# matrices created by create_n3d and hence to simplify the algorithm, but only the first N
# portion is used.
create_param3d(N::SVec3Int) =
    (s = (N.+1).data; (zeros(CFloat, s..., 3, 3), zeros(CFloat, s..., 3, 3)))  # 3 = numel(Axis)

# Below, zeros cannot be used instead of Array{T}, because zero for the type T may not be
# well-defined (e.g., T = Object3)
create_n3d(::Type{T}, N::SVec3Int) where {T} =
    (s = (N.+1).data; ((t->Array{T}(undef,t)).((s,s,s,s)), (t->Array{T}(undef,t)).((s,s,s,s))))  # Tuple24{Array{T,3}}


# Notes on ghost point transformation by boundary conditions:
#
# - Do not transform objects according to boundary conditions.  Periodization, for example,
# must be done by explicitly putting the translated object.  (In other words, the solver
# periodizes the whatever composition inside the domain, but composing the space inside the
# domain is the user's responsibility.)
#
# - Instead, transform ghost points back to the corresponding points inside the domain, and
# see which object is there.  Take that object to the ghost points.
#
# - Ghost points are usually copied from different non-ghost points by translation (for
# periodic BC) or reflection (for symmetry BC), but NOT if their transformed points are
# themselves (e.g., primal ghost points on PPC, and dual ghost points on PDC).  For the
# latter case, the fields of interest on the ghost points are still inferrable by other
# means than spatial transformation (i.e., the fields are zero at primal ghost points on PPC
# and dual ghost points on PDC).  Therefore, at the equation solving step we don't need to
# keep the degrees of freedom for those ghost points.  Still, to smooth material parameterts
# at non-ghost points around ghost boundaries, we need objects assigned to ghost points
# (because ghost points are the corners of the voxels centered at those non-ghost points).
# These ghost points are not copied from other non-ghost points by spatial transformation,
# so we need to assign objects to these ghost points.
#
# - When the ghost points are the copies of non-ghost points, both the fields and objects at
# the ghost points must be copied for consistency, not only one of them.  For example, if a
# periodic boundary condition is used, the primal grid point at the positive end must be
# copies of the primal grid point at the negative end, from which both the objects and
# fields must be copied to the ghost points.

# Notes on assignment:
#
# The goal is to "paint" grid points with a given object.  Initially I thought I would not
# need to paint ghost points independently, because I thought ghost points are by definition
# points whose properties are inferrable from other points.  This is not always the case.
# See the above "Notes on ghost point transformation by boundary conditions".
#
# If the ghost points are not exactly on symmetry BC, we can recover the ghost point arrays
# from somewhere else.  The algorithm below prepares grid point indices such that we can
# assign objects the same way regardless of the kinds of boundary conditions.

# Why do we need to set up so many arrays: param3d, obj3d, pind3d, oind3d?  In principle, it
# is sufficient to set up only param3d and obj3d, because in the smoothing algorithm we need
# to know pind and oind only inside a single voxel at a moment, and these two 2×2×2 arrays
# can be easily constructed from the obj3d array.
# However, retrieving pind and oind from obj3d point-by-point like this is very slow due to
# dynamic dispatch.  To reduce the amount of dynamic dispatch, oind3d and pind3d must be set
# up object-by-object rather than point-by-point.  This means that it is inevitable to
# construct pind3d and oind3d arrays.
function assign_param!(param3d::Tuple2{AbsArrComplex{5}},  # parameter array to set
                       obj3d::Tuple24{AbsArr{<:Object3,3}},  # object array to set
                       pind3d::Tuple24{AbsArr{ParamInd,3}},  # material parameter index array to set
                       oind3d::Tuple24{AbsArr{ObjInd,3}},  # object index array to set
                       ovec::AbsArr{<:Object3},  # object vector; later object overwrites earlier.
                       τl::Tuple23{AbsVecReal},  # field component locations transformed by boundary conditions
                       isbloch::SVec3Bool)
    # Circularly shift subscripts for Bloch boundary condition.  This makes sure τl[sub[w]]
    # is always sorted for all boundary conditions.  Sorted τl is necessary to use findfirst
    # and findlast in assign_param_obj!.
    #
    # Primal grid: ghost points are always at the positive end.
    # - Bloch: ghost points copy the values of non-ghost points at the negative end, so τl
    # must be circshifted to right in order to be sorted.
    # - symmetry: ghost points have own degrees of freedom at the positive end, so τl is
    # already sorted.
    # Therefore, the primal grid needs to be circshifted to right only for the Bloch BC.
    #
    # Dual grid: ghost points are always at the negative end.
    # - Bloch: ghost points copy the values of non-ghost points at the positive end, so τl
    # must be circshifted to left in order to be sorted.
    # - symmetry: ghost points copy the values of non-ghost points at the negative end, so τl is
    # already sorted.
    M = length.(τl[nPR])  # N+1
    sub = (map((m,b)->circshift(1:m, b), M, isbloch.data),
           map((m,b)->circshift(1:m, -b), M, isbloch.data))

    # Subscripts for param3d are a bit different.  param3d always have ghost cells at the
    # positive end, regardless of the material parameter type.  Therefore, if the above sub
    # puts ghost points at the positive end, the psub below does not have to be circshifted.
    # On the other hand, if the above sub puts ghost points at the negative end, the psub
    # below must be circshifted by +1 to push the ghost points to the negative end.
    #
    # Bottom line: prepare subscripts such that circshifted param3d and circshifted τl have
    # ghost points at the same indices.  Note that before circshift,
    # - param3d has ghost cells at the positive end,
    # - the primal grid has ghost points at the positive end, and
    # - the dual grid has ghost points at the negative end.
    # This means that before circshift, primal param3d and the primal τl have ghost points
    # at the same indices, whereas dual param3d and the dual τl don't.  Therefore, primal
    # param3d must be circshifted to right if the primal τl is circshifted to right.  (This
    # is to send param3d's ghost points existing at the positive end to the negative end).
    # Similarly, dual param3d must be circshifted to right if the dual τl is NOT circshifted
    # to left.  (This is also to send param3d's ghost points existing at the positive end to
    # the negative end, where the ghost points of the non-circshifted dual τl exist).
    psub = (map((m,b)->circshift(1:m, b), M, isbloch.data),  # primal grid
            map((m,b)->circshift(1:m, !b), M, isbloch.data))  # dual grid

    ## Perform assignment.
    for ngt = nPD
        gt = PD[ngt]
        ngt′ = alter(ngt)
        gt_cmp₀ = SVector(gt, gt, gt)
        for nw = 1:4
            # Set the grid types of the x-, y-, z-locations of Fw.
            gt_cmp = broadcast((k,w,g)->(k==w ? alter(g) : g), nXYZ, nw, gt_cmp₀)  # no change if nw = 4

            # Choose the circularly shifted subscripts to use.
            sub_cmp = t_ind(sub, gt_cmp)
            psub_cmp = t_ind(psub, gt_cmp)

            # Prepare the circularly shifted locations of the field components.
            τlcmp = view.(t_ind(τl,gt_cmp), sub_cmp)  # Tuple3{VecFloat}: locations of Fw = Uw or Vw

            # Prepare the circularly shifted viewes of various arrays to match the sorted
            # τl.  Even though all arrays are for same locations, param3d_cmp contains gt
            # material, whereas obj3d_cmp, pind3d_cmp, oind3d_cmp contain alter(gt)
            # material (so use ngt′ instead of ngt for them).
            param3d_cmp = view(param3d[ngt], psub_cmp..., nXYZ, nXYZ)
            obj3d_cmp = view(obj3d[ngt′][nw], sub_cmp...)
            pind3d_cmp = view(pind3d[ngt′][nw], sub_cmp...)
            oind3d_cmp = view(oind3d[ngt′][nw], sub_cmp...)

            # Set various arrays for the current component.
            assign_param_cmp!(gt, nw, param3d_cmp, obj3d_cmp, pind3d_cmp, oind3d_cmp, ovec, τlcmp)
        end
    end

    return nothing
end


# This is the innermost function that is customized for MaxwellFDM (i.e., it takes arrays
# of specific element types).
function assign_param_cmp!(gt::GridType,  # primal field (U) or dual field (V)
                           nw::Integer,  # w = X̂ (1), Ŷ (2), Ẑ (3), grid node (4)
                           param3d_cmp::AbsArrComplex{5},  # parameter array to set
                           obj3d_cmp::AbsArr{Object3,3},  # object array to set
                           pind3d_cmp::AbsArr{ParamInd,3},  # material parameter index array to set
                           oind3d_cmp::AbsArr{ObjInd,3},  # object index array to set
                           ovec::AbsVec{<:Object3},  # object vector; later object overwrites earlier.
                           τlcmp::Tuple3{AbsVecReal})  # location of field components
    for o = ovec  # last in ovec is last object added; see object.jl/add!
        # Retrieve shape once here, so that it can be passed to the function barrier.
        shape = o.shape

        # Prepare values to set.
        gt′ = alter(gt)
        param = matparam(o,gt)
        pind′ = paramind(o,gt′)
        oind = objind(o)

        # Set the values to the arrays.
        #
        # In fact, assign_val_shape! was designed to perform assignment on multiple arrays
        # at once.  This is meant to be achieved by passing a tuple named `arrays` below.
        # However, this is an inhomogeneous tuple, and retrieving entries from a tuple uses
        # too many allocations.  See the problem reported at
        # - https://github.com/JuliaLang/julia/issues/19850
        # - https://discourse.julialang.org/t/broadcasting-setindex-over-a-tuple-of-arrays-with-splatted-indices-is-slow/9641
        #
        # So, use the more compact code (commented below) once this issue is resolved.  The
        # current code performs Base.in(pt, Shape) more than necessary, but
        # benchmark/smoothing.jl demonstrates that the performance does not degrade much.
        assign_val_shape!(pind3d_cmp, pind′, shape, τlcmp)
        assign_val_shape!(oind3d_cmp, oind, shape, τlcmp)
        assign_val_shape!(obj3d_cmp, o, shape, τlcmp)
        if nw == 4
            assign_val_shape!(param3d_cmp, param, shape, τlcmp)
        else  # nw = 1, 2, 3
            assign_val_shape!(@view(param3d_cmp[:,:,:,nw,nw]), param[nw,nw], shape, τlcmp)
        end
        # arrays = (pind3d_cmp, oind3d_cmp, obj3d_cmp)
        # vals = (pind′, oind, o)
        # if nw == 4
        #     assign_val_shape!((arrays..., param3d_cmp), (vals..., param), shape, τlcmp)
        # else  # nw = 1, 2, 3
        #     assign_val_shape!((arrays..., @view(param3d_cmp[:,:,:,nw,nw])), (vals..., param[nw,nw]), shape, τlcmp)
        # end
    end

    return nothing
end


# Do we need an another wrapper function that selects τlcmp from gt and nw?  The downside of
# such a function is that avf already contains nw information.
# Maybe, I could create a function that selects τlcmp from gt and nw.


# assign_val_shape!(array::AbsArr{T,3}, val::T, shape::Shape{3}, τlcmp::Tuple3{AbsVecReal}) where {T} =
#     assign_val_shape!((array,), (val,), shape, τlcmp)
#
# assign_val_shape!(array::AbsArr{T,5}, val::AbsMat{T}, shape::Shape{3}, τlcmp::Tuple3{AbsVecReal}) where {T} =
#     assign_val_shape!((array,), (val,), shape, τlcmp)
#
# # Given a shape, assign value at the points within the shape.
# #
# # This function is written to handle a variable-length list of (array, value, function)'s in
# # assign_param_cmp, but the usage is suppressed for performance.
# function assign_val_shape!(arrays::Tuple,
#                            vals::Tuple,
#                            shape::Shape{3},
#                            τlcmp::Tuple3{AbsVecReal})
#     # Set the location indices of object boundaries.
#     @assert all(issorted.(τlcmp))
#     bn, bp = bounds(shape)  # (SVec3, SVec3)
#     subn = map((l,b) -> (n = findfirst(l.≥b); n==nothing ? 1 : n), τlcmp, bn)  # SVec3Int
#     subp = map((l,b) -> (n = findlast(l.≤b); n==nothing ? length(l) : n), τlcmp, bp)  # SVec3Int
#     I, J, K = map((nᵢ,nₑ) -> nᵢ:nₑ, subn, subp)  # SVec3{UnitRange{Int}}
#
#     if shape isa Box{3,9} && (shape::Box{3,9}).p == SMatrix{3,3,Float}(LinearAlgebra.I)  # shape is Cartesian box
#         # @info "haha1"
#         assign_val_range!.(arrays, vals, ((I,J,K),))
#         # @info "haha2"
#     else  # shape is not Cartesian box
#         # @info "haha3"
#         for k = K, j = J, i = I  # z-, y-, x-indices
#             pt = t_ind(τlcmp, i, j, k)
#             if pt ∈ shape
#                 @time assign_val!.(arrays, vals, ((i,j,k),))
#             end  # if pt ∈ shape
#         end  # for k = ..., j = ..., i = ...
#         # @info "haha4"
#     end  # if shape isa ...
#
#     return nothing
# end

# Do we need an another wrapper function that selects τlcmp from gt and nw?  The downside of
# such a function is that avf already contains nw information.
# Maybe, I could create a function that selects τlcmp from gt and nw.


assign_val_shape!(array::AbsArr{T,3}, val::T, shape::Shape{3}, τlcmp::Tuple3{AbsVecReal}) where {T} =
    assign_val_shape_impl!(array, val, shape, τlcmp)

assign_val_shape!(array::AbsArr{T,5}, val::AbsMat{T}, shape::Shape{3}, τlcmp::Tuple3{AbsVecReal}) where {T} =
    assign_val_shape_impl!(array, val, shape, τlcmp)

# Given a shape, assign value at the points within the shape.
#
# This function is written to handle a variable-length list of (array, value, function)'s in
# assign_param_cmp, but the usage is suppressed for performance.
function assign_val_shape_impl!(array::AbsArr{T},
                           val::Union{T,AbsMat{T}},
                           shape::Shape{3},
                           τlcmp::Tuple3{AbsVecReal}) where {T}
    # Set the location indices of object boundaries.
    @assert all(issorted.(τlcmp))
    bn, bp = bounds(shape)  # (SVec3, SVec3)
    subn = map((l,b) -> (n = findfirst(l.≥b); n==nothing ? 1 : n), τlcmp, bn.data)  # NTuple{3,Int}
    subp = map((l,b) -> (n = findlast(l.≤b); n==nothing ? length(l) : n), τlcmp, bp.data)  # NTuple{3,Int}
    I, J, K = map((nᵢ,nₑ) -> nᵢ:nₑ, subn, subp)  # NTuple{3,UnitRange{Int}}


    if shape isa Box{3,9} && (shape::Box{3,9}).p == LinearAlgebra.I  # shape is Cartesian box
        # @info "haha1"
        assign_val_range!(array, val, I, J, K)
        # @info "haha2"
    else  # shape is not Cartesian box
        # @info "haha3"
        for k = K, j = J, i = I  # z-, y-, x-indices
            pt = t_ind(τlcmp, i, j, k)
            if pt ∈ shape
                # @info "pt = $pt, shape = $shape, typeof(array) = $(typeof(array))"
                # @btime assign_val2!($array, $val, ($i,$j,$k))
                assign_val!(array, val, i, j, k)
            end  # if pt ∈ shape
        end  # for k = ..., j = ..., i = ...
        # @info "haha4"
    end  # if shape isa ...

    return nothing
end

# Could be named Base.setindex!, but didn't want this to be exported, so named different.
function assign_val!(array::AbsArr{T,3}, scalar::T, i::Integer, j::Integer, k::Integer) where {T}
    # @info "haha7"
    @inbounds array[i,j,k] = scalar
    # @info "haha8"
    return nothing
end

function assign_val!(array::AbsArr{T,5}, tensor::AbsMat{T}, i::Integer, j::Integer, k::Integer) where {T}
    for nc = nXYZ, nr = next2(nc)  # column- and row-indices
        # @info "haha9"
        @inbounds array[i,j,k,nr,nc] = tensor[nr,nc]
        # @info "haha10"
    end
    return nothing
end

function assign_val_range!(array::AbsArr{T,3}, scalar::T, i::AbsVecInteger, j::AbsVecInteger, k::AbsVecInteger) where {T}
    # @info "haha5"
    @inbounds array[i,j,k] .= Ref(scalar)
    # @info "haha6"
    return nothing
end

function assign_val_range!(array::AbsArr{T,5}, tensor::AbsMat{T}, i::AbsVecInteger, j::AbsVecInteger, k::AbsVecInteger) where {T}
    for nc = nXYZ, nr = next2(nc)  # column- and row-indices
        # @info "haha9"
        @inbounds array[i,j,k,nr,nc] .= tensor[nr,nc]
        # @info "haha10"
    end
    return nothing
end
