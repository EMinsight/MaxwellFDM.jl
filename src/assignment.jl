# Below, obj3d and oind3d are different: obj3d stores references to Object{3} instances,
# whereas oind3d stores object indices.  The reason for this is to handle the case where
# objects touching the negative-side boundary and objects touching the positive-side
# boundary are joined to form a single object by being wrapped around by the periodic
# boundary condition.  In that case, we assign the same object index to those distinct
# objects.
#
# I could enhance the assignment performance by constructing oid3d (not oind3d), which
# stores unique object IDs rather than the reference itself.  This array is different from
# oind3d in that it distinguishes objects joined by periodic boundary condition.  Then, I
# can retrieve objects from a map from this IDs to objects.  I could use oid3d in smoothing
# as well, by constructing an 8-vector oid3d_vxl, like pind3d_vxl and oind3d_vxl.
#
# The rationale for this trick is that assigning an Int to VecInt is faster than assigning a
# concrete Object (like Box) to Vector{Object{3}}.  In my test,
#   setindex!(::VecInt, ::Int, <index>)
# was about twice as fast as
#   setindex!(::Vector{Object{3}}, ::Box, <index>), though both took only a few nanoseconds.
#
# Currently the performance of assignment and smoothing seems satisfactory, so I will not
# pursue this extra optimization.  If I need to assign many objects, the situation may
# change and I may need to implement this optimization.

export create_param_array, create_p_storage, assign_param!, assign_val_shape!

# About the order of indices of param3d
#
# param3d has 5 indices: the first three are positional indices (i,j,k) and the last two are
# material parameter tensor component indices (v,w).  This mean param3d[i,j,k,:,:] is the
# 3×3 material tensor.
#
# This is different from the conventional indexing scheme used for material parametetrs.
# For example, we usually write the xy-component of ε at a location (i,j,k) as ε_xy[i,j,k],
# where the tesnor component subscripts v = x and w = y comes earlier than the positional
# indices (i,j,k).  Similarly, the x-component of the E-field is usually written E_x[i,j,k],
# not E[i,j,k]_x.
#
# However, when assigning the material parameters in the code in this file, we often assign
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


# About the material tensors stored as param3d
#
# param3d[ft][i,j,k,:,:] is a 3×3 tensor.  Suppose we are dealing with electric materials.
# Then this is the ε tensor.  It is important to note that not all the entries of this
# tensor are measured at the same physical locations.
#
# For example, if Ex, Ey, Ez are all tangential to primal grid planes (like in the standard
# Yee's grid construction), hten εxx, εyy, εzz are defined at (i+1/2,j,k), (i,j+1/2,k),
# (i,j,k+1/2), respectively.  Also, six εvw's with v ≠ w are defined at (i,j,k).


# Creates param3d to be of size N.+1 in the (i,j,k)-dimensions.  It is created with size N+1
# to match the sizes of other matrices created by create_p_storage and hence to simplify the
# algorithm, but only the first N portion is used.
#
# For length(N) = 2, the output param_array is indexd as param_array[i,j,v,w], where (i,j)
# is the grid cell location, and v and w are the row and column indices of the ncmp×ncmp tensorial
# material parameters (like the ε tensor and μ tensor).
create_param_array(N::SInt, ncmp::Int=3) = (s = (N.+1).data; zeros(CFloat, s..., ncmp, ncmp))  # (i,j,...) element is ncmp×ncmp tensor

# Below, zeros cannot be used instead of Array{T}, because zero for the type T may not be
# well-defined (e.g., T = Object{3})
#
# The output p_array is indexed as n3d[i,j,k,w].  For example, if w = X̂ and we are dealing
# with electric materials, then n3d[i,j,k,w] indicates some properties related to Ex[i,j,k].
# Below, we have n3d = obj3d, pind3d, oind3d, and n3d[i,j,k,X̂] is used to store the electric
# material properties evaluated at the Hx[i,j,k] (not Ex[i,j,k]) point.  These electric
# material properties that are seemingly evaluated at wrong locations are still related to
# Ex[i,j,k], because they are used to determine how to smooth the electric properties at the
# Ex[i,j,k] point in smoothing.jl.
create_p_storage(::Type{T}, N::SInt, ncmp::Int=4) where {T} = (s = (N.+1).data; Array{T}(undef, s..., ncmp))  # (i,j,...) element is ncmp-vector


# Notes on ghost point transformation by boundary conditions:
#
# - Do not transform objects according to boundary conditions.  Periodization, for example,
# must be done by explicitly putting the translated object.  (In other words, the solver
# periodizes the whatever composition inside the domain, but composing the space inside the
# domain is the user's responsibility.)
#
# - Instead, transform ghost points back to the corresponding points inside the domain, and
# see which object is there.  Assign that object to the ghost points.
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


# Notes on assignments and ghost points
#
# The goal is to "paint" grid points with a given object.  Initially I thought I would not
# need to paint ghost points independently, because I thought ghost points are by definition
# points whose properties are inferrable from other points.  This is not always the case.
# See the above "Notes on ghost point transformation by boundary conditions".
#
# If the ghost points are not exactly on symmetry BC, we can recover the ghost point arrays
# from somewhere else.  The algorithm below prepares grid point indices such that we can
# assign objects the same way regardless of the kinds of boundary conditions.

# Notes on the reasons for many arrays to set up
#
# Why do we need to set up so many arrays: param3d, obj3d, pind3d, oind3d?  In principle, it
# is sufficient to set up only param3d and obj3d, because in the smoothing algorithm we need
# to know pind and oind only inside a single voxel at a moment, and these two 2×2×2 arrays
# can be easily constructed from the obj3d array.
# However, retrieving pind and oind from obj3d point-by-point like this is very slow due to
# dynamic dispatch.  To reduce the amount of dynamic dispatch, oind3d and pind3d must be set
# up object-by-object rather than point-by-point.  This means that it is inevitable to
# construct pind3d and oind3d arrays.

# Unlike smoothing.jl/smooth_param! and param.jl/param3d2mat, assign_param! has to handle
# both electric and magnetic material properties simultaneously; see the comments inside the
# function body towards the end of the function.
function assign_param!(param3d::Tuple2{AbsArrComplex{5}},  # (electric, magnetic) parameter arrays to set
                       obj3d::Tuple2{AbsArr{<:Object{3,3,3},4}},  # (electric, magnetic) object arrays to set
                       pind3d::Tuple2{AbsArr{ParamInd,4}},  # (electric, magnetic) material parameter index arrays to set
                       oind3d::Tuple2{AbsArr{ObjInd,4}},  # (electric, magnetic) object index arrays to set
                       boundft::SVector{3,FieldType},  # boundary field type
                       objvec::AbsArr{<:Object{3,3,3}},  # object vector; later object overwrites earlier.
                       τl::Tuple23{AbsVecReal},  # field component locations transformed by boundary conditions
                       isbloch::SBool{3})  # boundary conditions
    # `sub` is the subscripts for obj3d, pind3d, oind3d.
    #
    # Store circularly shifted subscripts in sub for Bloch boundary condition, such that
    # τl[sub[w]] is always sorted for all boundary conditions.  Sorted τl is necessary to
    # use findfirst and findlast in assign_param_obj!.
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
    M = length.(τl[nPR])  # N.+1
    sub = (map((m,b)->circshift(1:m,b), M, isbloch.data),  # Tuple3{VecInt}: primal grid
           map((m,b)->circshift(1:m,-b), M, isbloch.data))  # Tuple3{VecInt}: dual grid

    # `psub` is the subscripts for param3d.
    #
    # The subscripts for param3d are a bit different from the subscripts for other arrays.
    # param3d always have ghost cells at the positive end, regardless of the grid type.  (To
    # be more precise, only the first N portion of param3d is used to construct the Maxwell
    # operator, and this used portion exclude the ghost points.)  Therefore, if the above
    # sub has ghost points' subscripts at the positive end, the psub below does not have to
    # be circshifted.  On the other hand, if the above sub has ghost points's subscripts at
    # the negative end, the psub below must be circshifted by +1 to push the ghost points to
    # the negative end.
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
    psub = (map((m,b)->circshift(1:m,b), M, isbloch.data),  # Tuple3{VecInt}: primal grid
            map((m,b)->circshift(1:m,!b), M, isbloch.data))  # Tuple3{VecInt}: dual grid

    ## Perform assignment.
    for nft = nEH
        ft = EH[nft]
        nft′ = alter(nft)
        ft′ = EH[nft′]

        ptensorvec = matparam(objvec, ft)  # Kp×Kp×length(objvec): used to set up param3d
        pind′vec = paramind(objvec, ft′)  # used to set up pind3d
        oindvec = objind(objvec)  # used to set up oind3d

        gt_cmp₀ = ft2gt.(ft, boundft)  # grid type of corners of voxel whose edges at ft lines
        for nw = 1:4
            # Set the grid types of the x-, y-, z-locations of Fw.
            gt_cmp = broadcast((k,w,g)->(k==w ? alter(g) : g), nXYZ, nw, gt_cmp₀)  # grid type of Fw; no change for nw = 4

            # Choose the circularly shifted subscripts to use.
            sub_cmp = t_ind(sub, gt_cmp)  # Tuple3{VecInt}
            psub_cmp = t_ind(psub, gt_cmp)  # Tuple3{VecInt}

            # Prepare the circularly shifted locations of the field components.
            τlcmp = view.(t_ind(τl,gt_cmp), sub_cmp)  # Tuple3{VecFloat}: locations of Fw = Uw or Vw

            # Prepare the circularly shifted viewes of various arrays to match the sorted
            # τl.  Even though all arrays are for same locations (τlcmp), param_cmp
            # contains ft material, whereas obj_cmp, pind_cmp, oind_cmp contain
            # alter(ft) material (so use nft′ instead of nft for them).
            #
            # Explanation.  We set up obj3d, pind3d, oind3d (param3d excluded) such that
            # they contain the voxel corner information to use while determining whether
            # param3d at the voxel centers needs smoothing or not.  For example, to
            # determine whether param3d defined at an Ew point needs smoothing or not, we
            # take the voxel centered at the Ew point and examine pind at the eight voxel
            # corners.  If all the pind's at the eight voxel corners are the same, we assume
            # that the voxel is filled with a uniform material parameter identified by the
            # pind.  This means even thought these voxel corners are the H-field points, we
            # still need to examine the electric material properties at these voxel corners.
            #
            # Conversely, the Ew points, at which we are going to set param3d, are the voxel
            # corners of the voxel centered at the Hw points.  So, at the Ew points we
            # would want to set the magnetic material properties in obj3d, pind3d, oind3d.
            # This is why for the same τlcmp locations we set up the ft-entry of param3d but
            # alter(ft)-entries of obj3d, pind3d, and oind3d.
            #
            # Note that view() below is used to get a circularly shifted version of the
            # array, not a portion of.
            obj_cmp = view(obj3d[nft′], sub_cmp..., nw)  # circularly shifted nw-component of obj3d[nft′]
            pind_cmp = view(pind3d[nft′], sub_cmp..., nw)  # circularly shifted nw-component of pind3d[nft′]
            oind_cmp = view(oind3d[nft′], sub_cmp..., nw)  # circularly shifted nw-component of oind3d[nft′]

            # Set various arrays for the current component.
            if nw == 4
                # Set the off-diagonal entires of param3d.
                param_cmp = view(param3d[nft], psub_cmp..., nXYZ, nXYZ)  # circularly shifted param3d[nft]; ndims = 5
                assign_param_cmp!(param_cmp, pind_cmp, oind_cmp, obj_cmp, ptensorvec, pind′vec, oindvec, objvec, τlcmp)
            else  # nw = 1, 2, 3
                # Set the diagonal entries of param3d.
                param_cmp = view(param3d[nft], psub_cmp..., nw, nw)  # circularly shifted (nw,nw)-component of param3d[nft]; ndims = 3
                pscalarvec = @view(ptensorvec[nw,nw,:])
                assign_param_cmp!(param_cmp, pind_cmp, oind_cmp, obj_cmp, pscalarvec, pind′vec, oindvec, objvec, τlcmp)
            end
        end
    end

    return nothing
end


# This is the innermost function that is customized for MaxwellFDM (i.e., it takes arrays
# of specific element types).
#
# The goal of this function is to set all the arrays' entries, where the arrays are
# param_cmp, obj_cmp, pind_cmp, oind_cmp, at the same locations (τlcmp).  For each
# shape we iterate over grid points and test if each point is included in the shape.
# Because there are so many grid points, the total time taken for these tests is somewhat
# substantial, so we don't want to test the same point again.  Therefore, for a point that
# is tested being included in the shape, we want to set up all the arrays that need to be
# set up.  If the point is an Ew point, it is where the electric material properties ε_ww is
# set up in param3d, but it is also where the magnetic material parameter index is set up in
# pind3d.  This is why the ft′ entry of pind3d was chosen as pind_cmp in the enclosing
# function, where as the ft entry of param3d was chosen as param_cmp there.  This also
# explains why the ft′ component was chosen as pind′ whereas the ft component was chosen as
# param in the present function.

# For nw = 1, 2, 3.
function assign_param_cmp!(param_cmp::AbsArrComplex{K},  # output material parameter array (array of scalars)
                           pind_cmp::AbsArr{ParamInd,K},  # output material parameter index array
                           oind_cmp::AbsArr{ObjInd,K},  # output object index array
                           obj_cmp::AbsArr{<:Object{K,Ke,Km},K},  # output object array
                           paramvec::AbsVecComplex,  # input material parameter vector (vector of scalars); length = length(objvec)
                           pind′vec::AbsVec{ParamInd},  # input material parameter index vector; length = length(objvec)
                           oindvec::AbsVec{ObjInd},  # input object index vector; length = length(objvec)
                           objvec::AbsVec{<:Object{K,Ke,Km}},  # input object vector; later object overwrites earlier.
                           τlcmp::NTuple{K,AbsVecReal}) where {K,Ke,Km}  # location of field components
    for n = 1:length(objvec)  # last in objvec is last object added; see object.jl/add!
        # Prepare values to set.
        param = paramvec[n]  # used to set up param_cmp
        pind′ = pind′vec[n]  # used to set up pind_cmp
        oind = oindvec[n]  # used to set up oind_cmp
        o = objvec[n]  # used to set up obj_cmp

        # Retrieve shape once here, so that it can be passed to the function barrier.
        shape = o.shape

        # Set the values to the arrays.
        #
        # Below, each of the four arrays pind_cmp, oind_cmp, obj_cmp, param_cmp is
        # individually set by assign_val_shape!.  This means that for each point p and shape
        # sh, the same test p ϵ sh is repeated four times.  This is inefficient.
        # An alternative implementation is to create a function similar to assign_val_shape!
        # that takes all these four arrays, test p ϵ sh once, and set all the arrays
        # simultaneously.  This may seem more efficient, but this turns out to perform a lot
        # of dynamic dispatches of assign_val!, because four different realizations of
        # assign_val! depending on the output array type are called for every p ∈ sh in the
        # for loop.
        #
        # A huge number of dynamic dispatches can be avoided by setting up each array
        # individually as below.  Then, the output array type of assign_val! is determined
        # at the interface of assign_val_shape!, so the dynamic dispatch of assign_val!
        # occurs only once per assign_val_shape!.
        assign_val_shape!(pind_cmp, pind′, shape, τlcmp)
        assign_val_shape!(oind_cmp, oind, shape, τlcmp)
        assign_val_shape!(obj_cmp, o, shape, τlcmp)

        # Assign parameters at Yee's field points where the v-component E and H affects
        # the w=v components of D and B.  This will use assign_val(..., scalar, ...) for
        # setting diagonal entries of the material parameter tensor.
        assign_val_shape!(param_cmp, param, shape, τlcmp)
    end

    return nothing
end

# For nw = 4.
function assign_param_cmp!(param_cmp::AbsArrComplex{L},  # output material parameter array (array of tensors); L = K+2
                           pind_cmp::AbsArr{ParamInd,K},  # output material parameter index array
                           oind_cmp::AbsArr{ObjInd,K},  # output object index array
                           obj_cmp::AbsArr{<:Object{K,Ke,Km},K},  # output object array
                           paramvec::AbsArrComplex{3},  # input material parameter vector (vector of tensors; paramvec[:,:,n] is tensor); length = length(objvec)
                           pind′vec::AbsVec{ParamInd},  # input material parameter index vector; length = length(objvec)
                           oindvec::AbsVec{ObjInd},  # input object index vector; length = length(objvec)
                           objvec::AbsVec{<:Object{K,Ke,Km}},  # input object vector; later object overwrites earlier.
                           τlcmp::NTuple{K,AbsVecReal}) where {K,Ke,Km,L}  # location of field components
    for n = 1:length(objvec)  # last in objvec is last object added; see object.jl/add!
        # Prepare values to set.
        param = @view(paramvec[:,:,n])  # used to set up param_cmp
        pind′ = pind′vec[n]  # used to set up pind_cmp
        oind = oindvec[n]  # used to set up oind_cmp
        o = objvec[n]  # used to set up obj_cmp

        # Retrieve shape once here, so that it can be passed to the function barrier.
        shape = o.shape

        # Set the values to the arrays.
        #
        # Below, each of the four arrays pind_cmp, oind_cmp, obj_cmp, param_cmp is
        # individually set by assign_val_shape!.  This means that for each point p and shape
        # sh, the same test p ϵ sh is repeated four times.  This is inefficient.
        # An alternative implementation is to create a function similar to assign_val_shape!
        # that takes all these four arrays, test p ϵ sh once, and set all the arrays
        # simultaneously.  This may seem more efficient, but this turns out to perform a lot
        # of dynamic dispatches of assign_val!, because four different realizations of
        # assign_val! depending on the output array type are called for every p ∈ sh in the
        # for loop.
        #
        # A huge number of dynamic dispatches can be avoided by setting up each array
        # individually as below.  Then, the output array type of assign_val! is determined
        # at the interface of assign_val_shape!, so the dynamic dispatch of assign_val!
        # occurs only once per assign_val_shape!.
        assign_val_shape!(pind_cmp, pind′, shape, τlcmp)
        assign_val_shape!(oind_cmp, oind, shape, τlcmp)
        assign_val_shape!(obj_cmp, o, shape, τlcmp)

        # Assign parameters at voxel corners where the v-component E and H affects the
        # w≠v components of D and B.  This will use assign_val(..., tensor, ...) for
        # setting off-diagonal entries of the material parameter tensor.
        assign_val_shape!(param_cmp, param, shape, τlcmp)
    end

    return nothing
end


# The only role of these wrappers assign_val_shape! is to limit the types of array and val
# to either (AbsArr{T,3}, T) or (AbsArr{T,5}, AbsMat{T}).  There is no performance benefit.
assign_val_shape!(array::AbsArr{T,K}, val::T, shape::Shape{K}, τlcmp::NTuple{K,AbsVecReal}) where {T,K} =
    assign_val_shape_impl!(array, val, shape, τlcmp)

assign_val_shape!(array::AbsArr{T,L}, val::AbsMat{T}, shape::Shape{K}, τlcmp::NTuple{K,AbsVecReal}) where {T,K,L} =
    assign_val_shape_impl!(array, val, shape, τlcmp)

# Given a shape, assign value at the points within the shape.
function assign_val_shape_impl!(array::AbsArr{T},
                                val::Union{T,AbsMat{T}},
                                shape::Shape{K},
                                τlcmp::NTuple{K,AbsVecReal}) where {T,K}
    # Set the location indices of object boundaries.
    @assert all(issorted.(τlcmp))
    bn, bp = bounds(shape)  # (SVector{K}, SVector{K})

    # Below, if shape is beyond the positive end of the domain, then findlast(l.≤b) = length(l)
    # but findfirst(l.≥b) = nothing.  To prevent any assignment, we have to replace nothing
    # with length(l)+1.  Then, the index range becomes length(l)+1:length(l), over which no
    # iteration occurs.  Similar consideration for shape beyond the negative end of the
    # domain.
    subn = map((l,b) -> (n = findfirst(l.≥b); n==nothing ? length(l)+1 : n), τlcmp, bn.data)  # NTuple{K,Int}
    subp = map((l,b) -> (n = findlast(l.≤b); n==nothing ? 0 : n), τlcmp, bp.data)  # NTuple{K,Int}
    CI = CartesianIndices(map((nᵢ,nₑ) -> nᵢ:nₑ, subn, subp))  # CartesianIndices{K}

    # Below, I think using the type assertion in (shape::Box{K,K*K}).p achieves type
    # stability of p, because Box{K,L}.p is a type of SMatrix{N,N,Float64,L}.
    if shape isa Box{K,K*K} && (shape::Box{K,K*K}).p == LinearAlgebra.I  # shape is Cartesian box
        assign_val!(array, val, CI)
    else  # shape is not Cartesian box
        for ci = CI
            pt = t_ind(τlcmp, ci)
            if pt ∈ shape
                assign_val!(array, val, ci)
            end  # if pt ∈ shape
        end  # for ci = ...
    end  # if shape isa ...

    return nothing
end

# Could be named Base.setindex!, but didn't want this to be exported, so named different.
function assign_val!(array::AbsArr{T,K}, scalar::T, ci::CartesianIndex{K}) where {T,K}
    @inbounds array[ci] = scalar
    return nothing
end

function assign_val!(array::AbsArr{T,K}, scalar::T, CI::CartesianIndices{K}) where {T,K}
    @inbounds array[CI] .= Ref(scalar)
    return nothing
end

# Similar to assign_val!(..., scalar, ...), but set the off-diagonal entries of a given
# tensor.
function assign_val!(array::AbsArr{T,L}, tensor::AbsMat{T}, ci::CartesianIndex{K}) where {T,K,L}  # L = K+2
    Nr, Nc = size(tensor)

    # Set the entries below the main diagonal.
    for nc = 1:Nc, nr = nc+1:Nr  # column- and row-indices
        @inbounds array[ci,nr,nc] = tensor[nr,nc]
    end

    # Set the entries above the main diagonal.
    for nc = 2:Nc, nr = 1:nc-1  # column- and row-indices
        @inbounds array[ci,nr,nc] = tensor[nr,nc]
    end

    return nothing
end

function assign_val!(array::AbsArr{T,L}, tensor::AbsMat{T}, CI::CartesianIndices{K}) where {T,K,L}  # L = K+2
    Nr, Nc = size(tensor)

    # Set the entries below the main diagonal.
    for nc = 1:Nc, nr = nc+1:Nr  # column- and row-indices
        @inbounds array[CI,nr,nc] .= tensor[nr,nc]
    end

    # Set the entries above the main diagonal.
    for nc = 2:Nc, nr = 1:nc-1  # column- and row-indices
        @inbounds array[CI,nr,nc] .= tensor[nr,nc]
    end

    return nothing
end
