# I could enhance the assignment performance by constructing oid3d (not oind3d), which stores
# unique object IDs rather than the reference itself.  This array is different from oind3d
# in that it distinguishes objects repeated by periodic boundary condition.  Then, I can
# retrieve objects from a map from this IDs to objects.  I could use oid3d in smoothing as
# well, by constructing an 8-vector oid3d_vxl, like pind3d_vxl and oind3d_vxl.
#
# The rationale for this trick is that assigning an Int to Vector{Int} is faster than
# assigning a concrete Object (like Box) to Vector{Object3}.  In my test,
#   setindex!(::Vector{Int}, ::Int, <index>)
# was about twice as fast as
#   setindex!(::Vector{Object3}, ::Box, <index>), though both took only a few nanoseconds.
#
# Currently the performance of assignment and smoothing seems satisfactory, so I will not
# pursue this extra optimization.  If I need to assign many objects, the situation may
# change and I may need to implement this optimization.

export create_param3d, create_n3d, assign_param!, assign_val_shape!

# Creates param3d of size N+1.  It is created with size N+1 to match the sizes of other
# matrices created by create_n3d and hence to simplify the algorithm, but only the first N
# portion is used.
create_param3d(N::SVec3Int) =
    (s = (N+1).data; (Array{CFloat}(s..., 3, 3), Array{CFloat}(s..., 3, 3)))  # 3 = numel(Axis)

create_n3d(::Type{T}, N::SVec3Int) where {T} =
    (s = (N+1).data; (Array{T}.((s,s,s,s)), Array{T}.((s,s,s,s))))  # Tuple24{Array{T,3}}


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
function assign_param!(param3d::Tuple2{AbsArr{CFloat,5}},  # parameter array to set
                       obj3d::Tuple24{AbsArr{<:Object3,3}},  # object array to set
                       pind3d::Tuple24{AbsArr{ParamInd,3}},  # material parameter index array to set
                       oind3d::Tuple24{AbsArr{ObjInd,3}},  # object index array to set
                       ovec::AbsArr{<:Object3},  # object vector; later object overwrites earlier.
                       τl::Tuple23{AbsVecReal},  # field component locations transformed by boundary conditions
                       ebc::SVector{3,EBC})
    # Circularly shift subscripts for BLOCH boundary condition.  This makes sure τl[sub[w]]
    # is always sorted for all boundary conditions.  Sorted τl is necessary to use findfirst
    # and findlast in assign_param_obj!.
    #
    # Primal grid: ghost points are always at the positive end.
    # - Bloch: ghost points copy the values of non-ghost points at the negative end, so τl
    # must be circshifted to right in order to be sorted.
    # - PPC: ghost points have own degrees of freedom at the positive end, so τl is already
    # sorted.
    # - PDC: ghost points copy the values of non-ghost points at the positive end, so τl is
    # already sorted.
    # Therefore, the primal grid needs to be circshifted to right only for the Bloch BC.
    #
    # Dual grid: ghost points are always at the negative end.
    # - Bloch: ghost points copy the values of non-ghost points at the positive end, so τl
    # must be circshifted to left in order to be sorted.
    # - PPC: ghost points copy the values of non-ghost points at the negative end, so τl is
    # already sorted.
    # - PDC: ghost points have own degrees of freedom at the negative end, so τl is already
    # sorted.
    M = length.(τl[nPR])  # N+1
    sub = (map((m,e)->circshift(1:m, e==BLOCH), M, ebc.data),
           map((m,e)->circshift(1:m, -(e==BLOCH)), M, ebc.data))

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
    psub = (map((m,e)->circshift(1:m, e==BLOCH), M, ebc.data),  # primal grid
            map((m,e)->circshift(1:m, e≠BLOCH), M, ebc.data))  # dual grid

    ## Perform assignment.
    for ngt = nPD
        gt = PD[ngt]
        ngt′ = alter(ngt)
        for nw = 1:4
            # Set the grid types of the x-, y-, z-locations of Fw.
            gt_cmp = SVector(gt, gt, gt)
            gt_cmp = broadcast((k,w,g)->(k==w ? alter(g) : g), nXYZ, nw, gt_cmp)  # no change if nw = 4

            # Choose the circularly shifted subscripts to use.
            sub_cmp = t_ind(sub, gt_cmp)
            psub_cmp = t_ind(psub, gt_cmp)

            # Prepare the circularly shifted locations of the field components.
            τlcmp = view.(t_ind(τl,gt_cmp), sub_cmp)  # Tuple3{Vector{Float}}: locations of Fw = Uw or Vw

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
                           nw::Integer,  # w = XX (1), YY (2), ZZ (3), grid node (4)
                           param3d_cmp::AbsArr{CFloat,5},  # parameter array to set
                           obj3d_cmp::AbsArr{Object3,3},  # object array to set
                           pind3d_cmp::AbsArr{ParamInd,3},  # material parameter index array to set
                           oind3d_cmp::AbsArr{ObjInd,3},  # object index array to set
                           ovec::AbsVec{<:Object3},  # object vector; later object overwrites earlier.
                           τlcmp::Tuple3{AbsVecFloat})  # location of field components
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


assign_val_shape!(array::AbsArr{T,3}, val::T, shape::Shape{3}, τlcmp::Tuple3{AbsVecFloat}) where {T} =
    assign_val_shape!((array,), (val,), shape, τlcmp)

assign_val_shape!(array::AbsArr{T,5}, val::AbsMat{T}, shape::Shape{3}, τlcmp::Tuple3{AbsVecFloat}) where {T} =
    assign_val_shape!((array,), (val,), shape, τlcmp)

# Given a shape, assign value at the points within the shape.
# Get a variable-length list of (array, value, function)'s.
function assign_val_shape!(arrays::Tuple,
                           vals::Tuple,
                           shape::Shape{3},
                           τlcmp::Tuple3{AbsVecFloat})
    # Set the location indices of object boundaries.
    assert(all(issorted.(τlcmp)))
    bn, bp = bounds(shape)  # (SVector{3}, SVector{3})
    subn = map((l,b) -> (n = findfirst(l.≥b); n==0 ? 1 : n), τlcmp, bn)  # SVec3Int
    subp = map((l,b) -> (n = findlast(l.≤b); n==0 ? length(l) : n), τlcmp, bp)  # SVec3Int
    I, J, K = map((nᵢ,nₑ) -> nᵢ:nₑ, subn, subp)  # SVector{3,UnitRange{Int}}

    if shape isa Box{3,9} && (shape::Box{3,9}).p == @SMatrix(eye(3))  # shape is Cartesian box
        assign_val!.(arrays, vals, ((I,J,K),))
    else  # shape is not Cartesian box
        for k = K, j = J, i = I  # z-, y-, x-indices
            pt = t_ind(τlcmp, i, j, k)
            if pt ∈ shape
                assign_val!.(arrays, vals, ((i,j,k),))
            end  # if pt ∈ shape
        end  # for k = ..., j = ..., i = ...
    end  # if shape isa ...

    return nothing
end


# Could be named Base.setindex!, but didn't want this to be exported, so named different.
assign_val!(array::AbsArr{T,3}, scalar::T, subs::Tuple3{S}) where {T,S<:Union{Integer,AbsVecInteger}} =
    (array[subs...] = scalar; nothing)

function assign_val!(array::AbsArr{T,5}, tensor::AbsMat{T}, subs::Tuple3{S}) where {T,S<:Union{Integer,AbsVecInteger}}
    for nc = nXYZ, nr = next2(nc)  # column- and row-indices
        array[subs..., nr, nc] = tensor[nr,nc]
    end
    return nothing
end
