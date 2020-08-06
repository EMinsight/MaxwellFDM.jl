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

export create_param_array, create_oind_array, assign_param!, assign_val_shape!

# About the order of indices of paramKd
#
# paramKd has 5 indices for K = 3: the first three are positional indices (i,j,k) and the
# last two are material parameter tensor component indices (v,w).  This mean paramKd[i,j,k,:,:]
# is the 3×3 material tensor.
#
# This is different from the conventional indexing scheme used for material parametetrs.
# For example, we usually write the xy-component of ε at a location (i,j,k) as ε_xy[i,j,k],
# where the tesnor component subscripts v = x and w = y comes earlier than the positional
# indices (i,j,k).  Similarly, the x-component of the E-field is usually written E_x[i,j,k],
# not E[i,j,k]_x.
#
# However, when assigning the material parameters in the code in this file, we often assign
# the same value to a block of paramKd.  This block is chosen differently for different
# material tensor components, because different tensor components are evaluated at different
# locations in Yee's grid.  Therefore, in this assignment, we first fix the material tensor
# components v and w, determine the range of (i,j,k) to which the same material parameter
# value to assign, and perform the assignment.
#
# Now, if we index paramKd as paramKd[v,w,i,j,k], v and w are the fastest-varing indices.
# Therefore, fixing them and assigning for a contiguous range of (i,j,k) does not actually
# assign in a contiguous memory block.  On the other hand, if we index paramKd as
# paramKd[i,j,k,v,w], for fixed v and w a contiguous range of (i,j,k) actually assigns in
# a more contiguous memory block.  (The entries with adjacent i for the same j,k,v,w are
# actually adjacent in the memory space.)  Therefore, the latter indexing scheme results in
# faster performance.  In my experiment of assigning 8 objects in 3D, the latter turns out
# to be about 30% faster.
#
# For this reason, I use the less conventional indexing scheme of paramKd[i,j,k,v,w].


# About the material tensors stored as paramKd
#
# paramKd[i,j,k,:,:] is a 3×3 tensor.  Suppose we are dealing with electric materials. Then
# this is the ε tensor.  It is important to note that not all the entries of this tensor are
# measured at the same physical locations.
#
# For example, if Ex, Ey, Ez are all tangential to primal grid planes (like in the standard
# Yee's grid construction), then εxx, εyy, εzz are defined at (i+½,j,k), (i,j+½,k), (i,j,k+½),
# respectively.  Also, six εvw's with v ≠ w are defined at (i,j,k).


# Create paramKd to be of size N.+1 in the (i,j,k)-dimensions.  In order to simplify the
# assignment algorithm, paramKd is created with size N+1 to match the sizes of the arrays
# created by create_oind_arrays.
#
# For length(N) = 2, the output paramKd is indexd as paramKd[i,j,v,w], where (i,j) is the
# grid cell location, and v and w are the row and column indices of the ncmp×ncmp tensorial
# material parameters (like the ε tensor and μ tensor).
create_param_array(N::SInt, ncmp::Int=3) = (s = (N.+1).data; zeros(CFloat, s..., ncmp, ncmp))  # (i,j,...) element is ncmp×ncmp tensor

#  Create an array to store object indices.
create_oind_array(N::SInt) = zeros(ObjInd, (N.+1).data)

# Notes on ghost point transformation by boundary conditions
#
# - First of all, do not transform objects according to boundary conditions.  Periodization,
# for example, must be done by explicitly putting the translated object.  (In other words,
# the solver periodizes the whatever composition inside the domain, but composing the space
# inside the domain is the user's responsibility.)
#
# - Transform ghost points back to the corresponding points inside the domain, and see which
# object is there.  Assign that object to the ghost points.
#
# - Ghost points are usually copied from different non-ghost points by translation (for
# the periodic BC) or reflection (for the symmetry BC), but NOT if their transformed points
# are themselves (e.g., primal ghost points on the symmetry BC).  For the latter case, the
# fields of interest on the ghost points are still inferrable by other means than spatial
# transformation (i.e., the fields are zero at primal ghost points on the symmetry BC).
# Therefore, at the equation solving step we don't need to keep the degrees of freedom for
# those ghost points.  Still, to smooth material parameterts at non-ghost points around
# ghost boundaries, we need objects assigned to ghost points (because ghost points are the
# corners of the voxels centered at those non-ghost points).  These ghost points are not
# copied from other non-ghost points by spatial transformation, so we need to assign objects
# to these ghost points.
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
# If the ghost points are not exactly on the symmetry BC, we can recover the ghost point
# arrays from somewhere else.  The algorithm below prepares grid point indices such that we can
# assign objects the same way regardless of the kinds of boundary conditions.

# Notes on setting object indices instead of Objects themselves
#
# Why do we set up an array oindKd′ of object indices rather than an array of Objects
# themselves?  By setting up an array of object indices, we have to pass the maps from
# object indices to shapes (oind2shp), to parameter indices (oind2pind), and also a map from
# parameter indices to material parameters (pind2matprm) in order to get the shape and
# material parameter at each grid point.  If we set up an array of Objects, we would not
# need to pass all these maps, because an Object already contains the shape and material
# parameter.  Still, there are a few reasons to pass an array of object indices instead of
# an array of objects:
#
# - Dispatching the shape and material paramter from an Object turns out to be slow, because
# Objects have different concrete types depending on the type of the Shape it contains.
# This makes the multiple dispatch to kick in, which is slow.
#
# - Setting the argumenty types Base types (e.g., AbstractArray) rather custom types (e.g.,
# Object in the present context) increases the reusability of the code.  For example, there
# are cases where we want to pass a submatrix of the 3×3 material parameter tensor (rather
# than the material parameter tensor itself), because we want to treat the transverse and
# longitudinal components of the material parameter tensor separately (e.g., in the
# waveguide mode solver).  By taking the elementary types, this can be achieved without
# changing the assign_param! function, because we can create a vector of 2×2 submatrices of
# the material parameter tensors outside assign_param! and pass it.  If assign_param! took
# an array of Objects, handling such cases would be difficult because we created an Object
# with a 3×3 material parameter tensor in this case.

# Unlike smoothing.jl/smooth_param! and param.jl/param3d2mat, the paramKd and oindKd′ passed
# to assign_param! are usually (but not always; see benchmark/smoothing2d.jl) for different
# field types (e.g., paramKd for E-field and oindKd′ for the H-field).
#
# Explanation.  We set up oindKd′ such that it contains the voxel corner information to use
# in smoothing for determining whether paramKd at the voxel centers needs smoothing or not.
# For example, to determine whether paramKd for K = 3 defined at an Ew point needs smoothing
# or not, we take the voxel centered at the Ew point and examine the material parameter
# indices (pind) at the eight voxel corners.  If all the pind's at the eight voxel corners
# are the same, we assume that the voxel is filled with a uniform material parameter
# identified by the pind.  This means even though these voxel corners are the H-field points,
# we still need to examine the electric material properties at these voxel corners.
#
# Conversely, the Ew points, at which we are going to set paramKd, are the voxel corners of
# the voxel centered at the Hw points.  So, at the Ew points we want to set oindKd′ to use
# for smoothing the magnetic material parameterts.  This is why we need to pass paramKd and
# oindKd′ for different field types.  (Note that we assign paramKd and oindKd′ at the same
# locations, because we set them up point-by-point.)

# About Kp₊₁⏐₁
#
# The number of entries in the tuple oindKd′ is Kp₊₁⏐₁, which is Kp+1 or 1.  Note that Kp+1
# cannot be 1 because Kp ≥ 1.  Therefore, by checking Kp₊₁⏐₁ = 1 or not, we can check if it
# is 1 or Kp+1.
#
# Kp₊₁⏐₁ is neither K or Kp.  It is the parameter to control the assignment behavior.  If
# the different components of the material parameter tensor need to be evaluated at
# different locations within a single Yee cell, the number of locations is Kp+1 in all cases
# I can think of, where the Kp locations are where the Kp diagonal entries are evaluated,
# and the additional 1 location is where all the off-diagonal entries are evaluated.  This
# covers many standard cases, such as ε and μ of the standard 3D problems, and also ε of the
# 2D TE and μ of the 2D TM problems.

# Below, primed variables indicate that the variables may not be of the same material type
# as the non-primed variables.  Compare this with the usage of the prime in smoothing.jl
# that indicates voxel corner quantities instead of voxel center quantities.
function assign_param!(paramKd::AbsArrComplex{K₊₂},  # electric (magnetic) parameter array to set; K₊₂ = K+2, where 2 is rank of material tensor
                       oindKd′::NTuple{Kp₊₁⏐₁,AbsArr{ObjInd,K}},  # object index arrays to set at electric (magnetic) parameter locations; Kp₊₁⏐₁ = Kp+1 or 1
                       gt₀::SVector{K,GridType},  # grid type of voxel corners; generated by ft2gt.(ft, boundft)
                       oind2shp::AbsVec{Shape{K,K²}},  # map from oind to shape; K² = K^2
                       oind2pind::AbsVec{ParamInd},  #  map from oind to electric (magnetic) material parameter index
                       pind2matprm::AbsVec{SSComplex{Kp,Kp²}},  # map from pind to electric (magnetic) material parameters; Kp² = Kp^2
                       τl::Tuple2{NTuple{K,AbsVecReal}},  # field component locations transformed by boundary conditions
                       isbloch::SBool{K}  # boundary conditions
                      ) where {K,Kp,K²,Kp²,K₊₂,Kp₊₁⏐₁}
    @assert(K²==K^2 && Kp²==Kp^2 && K₊₂==K+2 && (Kp₊₁⏐₁==Kp+1 || Kp₊₁⏐₁==1))
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
    @inbounds M = length.(τl[nPR])  # N.+1
    sub = (map((m,b)->circshift(1:m,b), M, isbloch.data),  # NTuple{K,VecInt}: primal grid
           map((m,b)->circshift(1:m,-b), M, isbloch.data))  # NTuple{K,VecInt}: dual grid

    # `psub` is the subscripts for paramKd.
    #
    # The subscripts for paramKd are a bit different from the subscripts for other arrays.
    # paramKd always have ghost cells at the positive end, regardless of the grid type.  (To
    # be more precise, only the first N portion of paramKd is used to construct the Maxwell
    # operator, and this used portion exclude the ghost points.)  Therefore, if the above
    # sub has ghost points' subscripts at the positive end, the psub below does not have to
    # be circshifted.  On the other hand, if the above sub has ghost points's subscripts at
    # the negative end, the psub below must be circshifted by +1 to push the ghost points to
    # the negative end.
    #
    # Bottom line: prepare subscripts such that circshifted paramKd and circshifted τl have
    # ghost points at the same indices.  Note that before circshift,
    # - paramKd has ghost cells at the positive end,
    # - the primal grid has ghost points at the positive end, and
    # - the dual grid has ghost points at the negative end.
    # This means that before circshift, primal paramKd and the primal τl have ghost points
    # at the same indices, whereas dual paramKd and the dual τl don't.  Therefore, primal
    # paramKd must be circshifted to right if the primal τl is circshifted to right.  (This
    # is to send paramKd's ghost points existing at the positive end to the negative end).
    # Similarly, dual paramKd must be circshifted to right if the dual τl is NOT circshifted
    # to left.  (This is also to send paramKd's ghost points existing at the positive end to
    # the negative end, where the ghost points of the non-circshifted dual τl exist).
    psub = (map((m,b)->circshift(1:m,b), M, isbloch.data),  # Tuple3{VecInt}: primal grid
            map((m,b)->circshift(1:m,!b), M, isbloch.data))  # Tuple3{VecInt}: dual grid

    ## Perform assignment for Kp₊₁⏐₁ == 1.
    if Kp₊₁⏐₁ == 1
        sub_cmp = t_ind(sub, gt₀)  # Tuple3{VecInt}
        psub_cmp = t_ind(psub, gt₀)  # Tuple3{VecInt}
        τlcmp = view.(t_ind(τl,gt₀), sub_cmp)  # Tuple3{VecFloat}: locations of Fw = Uw or Vw

        oind′_cmp = view(oindKd′[1], sub_cmp...)  # circularly shifted oindKd′; note oindKd′ has only one entry

        # Set oindKd′ and on- and off-diagonal entires of paramKd.
        param_cmp = view(paramKd, psub_cmp..., 1:Kp, 1:Kp)  # circularly shifted paramKd; ndims = K+2
        assign_param_cmp!(param_cmp, oind′_cmp, oind2shp, oind2pind, pind2matprm, τlcmp)  # this does nothing for Kp = 1

        return nothing
    end

    ## Perform assignment for Kp₊₁⏐₁ == Kp+1.
    @assert(Kp₊₁⏐₁ == Kp+1)
    for nw = Kp+1:-1:1  # nw = Kp+1 is handled first to allow remaining nw's to partially overwrite
        # Set the grid types of the x-, y-, z-locations of Fw.
        gt_cmp = gt_w(nw, gt₀)

        # Choose the circularly shifted subscripts to use.
        sub_cmp = t_ind(sub, gt_cmp)  # Tuple3{VecInt}
        psub_cmp = t_ind(psub, gt_cmp)  # Tuple3{VecInt}

        # Prepare the circularly shifted locations of the field components.
        τlcmp = view.(t_ind(τl,gt_cmp), sub_cmp)  # Tuple3{VecFloat}: locations of Fw = Uw or Vw

        # Note that view() below is used to get a circularly shifted version of the
        # array, not a portion of.
        oind′_cmp = view(oindKd′[nw], sub_cmp...)  # circularly shifted nw-component of oindKd′

        if nw == Kp+1  # w = Ô
            # Set oindKd′ and the on- and off-diagonal entires of paramKd.
            param_cmp = view(paramKd, psub_cmp..., 1:Kp, 1:Kp)  # circularly shifted paramKd; ndims = K+2
            assign_param_cmp!(param_cmp, oind′_cmp, oind2shp, oind2pind, pind2matprm, τlcmp)  # this does nothing for Kp = 1
        else  # w = X̂, Ŷ, Ẑ for K = 3
            # Overwrite the diagonal entries of paramKd.
            param_cmp = view(paramKd, psub_cmp..., nw, nw)  # circularly shifted (nw,nw)-component of paramKd; ndims = 3
            pind2matprm_w = [mp[nw,nw] for mp = pind2matprm]
            assign_param_cmp!(param_cmp, oind′_cmp, oind2shp, oind2pind, pind2matprm_w, τlcmp)
        end
    end

    return nothing
end


# This is the innermost function that is customized for MaxwellFDM (i.e., it takes arrays
# of specific element types).
#
# The goal of this function is to set all the arrays' entries, where the arrays are
# param_cmp and oind′_cmp at the same locations (τlcmp).  For each shape we iterate over
# grid points and test if each point is included in the shape.  Because there are so many
# grid points, the total time taken for these tests is somewhat substantial, so we don't
# want to test the inclusion of the same point again.  Therefore, for a point that is
# confirmed included in the shape, we want to set up all the arrays that can be set up.
# If the point is an Ew point, it is where the electric material properties ε_ww is set up
# in paramKd, but it is also where the object index for smoothing magnetic material
# parameters is set up in oindKd′, because this Ew point is the corners of the voxel whose
# center is at an Hw point, and we use the material parameter index at these voxel corners
# to determine whether the voxel centered at the Hw point is filled with a single magnetic
# material.
#
# For w = X̂, Ŷ, Ẑ
assign_param_cmp!(param_cmp::AbsArrComplex{K},  # output material parameter array (array of scalars)
                  oind′_cmp::AbsArr{ObjInd,K},  # output object index array
                  oind2shp::AbsVec{Shape{K,K²}},  # input map from oind to shape
                  oind2pind::AbsVec{ParamInd},  # input map from oind to pind
                  pind2matprm::AbsVecComplex,  #  input map from pind to material parameter
                  τlcmp::NTuple{K,AbsVecReal}) where {K,K²} = # location of field components
    assign_param_cmp_impl!(param_cmp, oind′_cmp, oind2shp, oind2pind, pind2matprm, τlcmp)

# For w = Ô.
# Below, note that the element type of pint2matprm, SSComplex{Kp,Lp}, is always a matrix
# (i.e., rank-2 tensor), regardless of the value of Kp.  For example, if we are assigning μz
# in a 2D TE problem, μz is a type of SSComplex{1,1}.  This is why M = K + 2, where 2 is the
# rank of the material parameter tensor, regardless of the value of Kp.
assign_param_cmp!(param_cmp::AbsArrComplex{K₊₂},  # output material parameter array (array of tensors); M = K+2
                  oind′_cmp::AbsArr{ObjInd,K},  # output object index array
                  oind2shp::AbsVec{Shape{K,K²}},  # input map from oind to shape
                  oind2pind::AbsVec{ParamInd},  # input map from oind to pind
                  pind2matprm::AbsVec{SSComplex{Kp,Kp²}},  #  input map from pind to material parameter
                  τlcmp::NTuple{K,AbsVecReal}) where {K,Kp,K²,Kp²,K₊₂} = # location of field components
    assign_param_cmp_impl!(param_cmp, oind′_cmp, oind2shp, oind2pind, pind2matprm, τlcmp)

function assign_param_cmp_impl!(param_cmp, oind′_cmp, oind2shp, oind2pind, pind2matprm, τlcmp)
    for oind = 1:length(oind2shp)  # last index of oind2shp is index of last object; see object.jl/add!
        # Prepare values to set.
        @inbounds shp = oind2shp[oind]
        @inbounds pind = oind2pind[oind]
        @inbounds param = pind2matprm[pind]

        # Set the values to the arrays.
        #
        # Below, each of the two arrays oind′_cmp and param_cmp is individually set by
        # assign_val_shape!.  This means that for each point p and shape sh, the same test
        # p ϵ sh is repeated two times.  This is inefficient. An alternative implementation
        # is to create a function similar to assign_val_shape! that takes all these two
        # arrays, test p ϵ sh once, and set all the arrays simultaneously.  This may seem
        # more efficient, but this turns out to perform a lot of dynamic dispatches of
        # assign_val!, because two different realizations of assign_val! depending on the
        # output array type are called for every p ∈ sh in the for loop.
        #
        # A huge number of dynamic dispatches can be avoided by setting up each array
        # individually as below.  Then, the output array type of assign_val! is determined
        # at the interface of assign_val_shape!, so the dynamic dispatch of assign_val!
        # occurs only once per assign_val_shape!.
        assign_val_shape!(oind′_cmp, ObjInd(oind), shp, τlcmp)

        # Assign parameters at Yee's field points where the v-component E and H affects
        # the w=v components of D and B.  This will use assign_val(..., scalar, ...) for
        # setting diagonal entries of the material parameter tensor.
        assign_val_shape!(param_cmp, param, shp, τlcmp)
    end

    return nothing
end


# The only role of these wrappers assign_val_shape! is to limit the types of array and val
# to either (AbsArr{T,3}, T) or (AbsArr{T,5}, AbsMat{T}).  There is no performance benefit.
assign_val_shape!(array::AbsArr{T,K}, val::T, shape::Shape{K}, τlcmp::NTuple{K,AbsVecReal}) where {T,K} =
    assign_val_shape_impl!(array, val, shape, τlcmp)

assign_val_shape!(array::AbsArr{T,K₊₂}, val::SMatrix{Kp,Kp,T}, shape::Shape{K}, τlcmp::NTuple{K,AbsVecReal}) where {T,K,Kp,K₊₂} =
    assign_val_shape_impl!(array, val, shape, τlcmp)

# Given a shape, assign value at the points within the shape.
function assign_val_shape_impl!(array::AbsArr{T},
                                val::Union{T,SMatrix{Kp,Kp,T}},
                                shape::Shape{K},
                                τlcmp::NTuple{K,AbsVecReal}) where {T,K,Kp}
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

# Similar to assign_val!(..., scalar, ...), but set the on- and off-diagonal entries of a given
# tensor.
function assign_val!(array::AbsArr{T,K₊₂}, tensor::SMatrix{Kp,Kp,T}, ci::CartesianIndex{K}) where {T,K,Kp,K₊₂}  # K₊₂ = K+2
    for c = 1:Kp, r = 1:Kp  # column- and row-indices
        @inbounds array[ci,r,c] = tensor[r,c]
    end

    return nothing
end

function assign_val!(array::AbsArr{T,K₊₂}, tensor::SMatrix{Kp,Kp,T}, CI::CartesianIndices{K}) where {T,K,Kp,K₊₂}  # K₊₂ = K+2
    for c = 1:Kp, r = 1:Kp, ci = CI  # column- and row-indices and location index
        @inbounds array[ci,r,c] = tensor[r,c]
    end

    return nothing
end
