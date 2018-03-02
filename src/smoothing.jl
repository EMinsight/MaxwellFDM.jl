# How do we handle TF/SF?  Need a capability to do subpixel smoothing only inside some region.

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

export create_param3d, create_n3d, assign_param!, smooth_param!

# See http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array.
@inline function swap8!(ind::AbsVecInteger, v::AbsVec, i::Integer, j::Integer)
    @inbounds if v[ind[j]] < v[ind[i]]
        @inbounds ind[i], ind[j] = ind[j], ind[i]
    end
end
@inline function sort8!(ind::AbsVecInteger, v::AbsVec)
    # Initialize ind.
    @simd for n = 1:8
        @inbounds ind[n] = n
    end

    # Sort ind.
    swap8!(ind, v, 1, 2); swap8!(ind, v, 3, 4); swap8!(ind, v, 1, 3); swap8!(ind, v, 2, 4);
    swap8!(ind, v, 2, 3); swap8!(ind, v, 5, 6); swap8!(ind, v, 7, 8); swap8!(ind, v, 5, 7);
    swap8!(ind, v, 6, 8); swap8!(ind, v, 6, 7); swap8!(ind, v, 1, 5); swap8!(ind, v, 2, 6);
    swap8!(ind, v, 2, 5); swap8!(ind, v, 3, 7); swap8!(ind, v, 4, 8); swap8!(ind, v, 4, 7);
    swap8!(ind, v, 3, 5); swap8!(ind, v, 4, 6); swap8!(ind, v, 4, 5)

    return ind
end
@inline function countdiff(ind::AbsVecInteger, v::AbsVec)
    c, nₑ = 1, 9
    @simd for n = 2:8
        @inbounds v[ind[n-1]] ≠ v[ind[n]] && (c += 1; nₑ = n)
    end
    return c, nₑ  # nₑ: last n where v[ind[n]] changed
end

const NOUT_VXL =  # vectors from corners to center
[SVector(1.,1.,1.), SVector(-1.,1.,1.), SVector(1.,-1.,1.), SVector(-1.,-1.,1.),
SVector(1.,1.,-1.), SVector(-1.,1.,-1.), SVector(1.,-1.,-1.), SVector(-1.,-1.,-1.)]

# This function creates param3d of size N+1, but note that only the first N portion is used.
# It is created with size N+1 to match the sizes of other matrices and hence to simplify the
# algorithm.
create_param3d(N::SVec3Int) =
    (s = (N+1).data; (Array{CFloat}(s..., 3, 3), Array{CFloat}(s..., 3, 3)))  # 3 = numel(Axis)

create_n3d(::Type{T}, N::SVec3Int) where {T} =
    (s = (N+1).data; (Array{T}.((s,s,s,s)), Array{T}.((s,s,s,s))))  # Tuple24{Array{T,3}}

# Overall smoothing algorithm:
# - Assign obj, pind, oind to arrays object-by-object.
#     - For the locations of grid points to assign the objects to, use lcmp created considering BC (τlcmp).
# - Using pind and oind, determine voxels to perform subpixel smoothing.
# - Inside each of the voxels, figure out the foreground object with which subpixel smoothing is performed.
#     - Iterating voxel corners, find the object that is put latest.  That object is the foreground object.
#     - Complication occurs when the voxel corner assigned with the foreground object is outside the domain boundary.
#         - In that case, to calculate nout and rvol properly, we have to move the voxel center.  (This is effectively the same as moving the foreground object.)
#             - If the corner is outside the periodic boundary, translate the voxel center before surfpt_nearby (using ∆fg).
#             - If the corner is outside the symmetry boundary, zero the nout componnent normal to the boundary (using σvxl).
#             - ∆fg and σvxl can be obtained simply by looking at the indices and BC.

# Notes on ghost point transformation by boundary conditions:
#
# - Do not transform objects according to boundary conditions.  Periodization, for example,
# must be done by explicitly putting the translated object.  (In other words, the solver
# periodizes the whatever composition inside the domain, but composing the space inside the
# domain is the user's responsibility.)
#
# - Instead, transform ghost points back to the corresponding points inside the domain, and
# see which object is there.  Take that object to point the ghost points.
#
# - Ghost points are usually copied from different non-ghost points by translation (for
# periodic BC) or reflection (for symmetry BC), but they are not if their transformed points
# are themselves (e.g., primal ghost points on PPC, and dual ghost points on PDC).  For the
# latter case, the fields of interest are still inferrable by means other than spatial
# transformation (i.e., the field are zero at primal ghost points on PPC and dual ghost
# points on PDC), so in the equation solving step we don't have to keep degrees of freedom
# for these ghost points.  However, to smooth material parameterts at non-ghost points
# around ghost boundaries, we need objects assigned to ghost points.  If these ghost points
# are not copied from other non-ghost points by spatial transformation, we need to assign
# objects to these ghost points.
#
# - When the ghost points are copies of non-ghost points, both the fields and objects at the
# ghost must be copied for consistency, not only one of them.  For example, if a periodic
# boundary condition is used, the primal grid point at the positive end must be copies of
# the primal grid point at the negative end, from which both the objects and fields at the
# ghost points must be copied.

# Notes on assignment:
#
# The goal is to "paint" grid points with a given object.  Initially I thought I would not
# need to paint ghost points independently, because I thought ghost points are by definition
# points whose properties are inferrable from other points.  This is not always the case.
# See the above "Notes on ghost point transformation by boundary conditions".
#
# If the ghost points are not exactly on the symmetry boundary condition, we can recover the
# ghost point arrays from somewhere else.  The algorithm below prepares grid point indices
# such that we can assign objects the same way regardless of the kinds of boundary conditions.

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
    M = length.(τl[nPR])  # N+1
    sub = (map((m,e)->circshift(1:m, e==BLOCH), M, ebc.data),
           map((m,e)->circshift(1:m, -(e==BLOCH)), M, ebc.data))

    # Subscripts for param3d are a bit different.  param3d always have ghost cells at the
    # positive end, regardless of the parameter type.  Therefore, if the above sub puts
    # ghost points at the positive end, the psub below does not have to be circshifted.  On
    # the other hand, if the above sub puts ghost points at the negative end, the psub below
    # must be circshifted by +1.
    #
    # Bottom line: prepare subscripts such that param3d and other arrays have ghost cells
    # at the same subscripts.  Note that before circshift,
    # - param3d has ghost cells at the positive end,
    # - the primal grid has ghost points at the positive end, and
    # - the dual grid has ghost points at the negative end.
    # This means that before circshift, param3d and the primal grid have the same ghost cell
    # location, whereas param3d and the dual grid don't.  Therefore, param3d must be
    # circshifted if the primal grid is circshifted or the dual grid is NOT circshifted.
    psub = (map((m,e)->circshift(1:m, e==BLOCH), M, ebc.data),
            map((m,e)->circshift(1:m, e≠BLOCH), M, ebc.data))

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
            param3d_cmp = view(param3d[ngt], psub_cmp..., 1:3, 1:3)
            obj3d_cmp = view(obj3d[ngt′][nw], sub_cmp...)
            pind3d_cmp = view(pind3d[ngt′][nw], sub_cmp...)
            oind3d_cmp = view(oind3d[ngt′][nw], sub_cmp...)

            # Set various arrays for the current component.
            assign_param_cmp!(gt, nw, param3d_cmp, obj3d_cmp, pind3d_cmp, oind3d_cmp, ovec, τlcmp)
        end
    end

    return nothing
end

function assign_param_cmp!(gt::GridType,  # primal field (U) or dual field (V)
                           nw::Int,  # w = XX (1), YY (2), ZZ (3), grid node (4)
                           param3d_cmp::AbsArr{CFloat,5},  # parameter array to set
                           obj3d_cmp::AbsArr{Object3,3},  # object array to set
                           pind3d_cmp::AbsArr{ParamInd,3},  # material parameter index array to set
                           oind3d_cmp::AbsArr{ObjInd,3},  # object index array to set
                           ovec::AbsVec{<:Object3},  # object vector; later object overwrites earlier.
                           τlcmp::Tuple3{AbsVec{Float}})  # location of field components
    for n = length(ovec):-1:1  # last in ovec is first object added; see shape.jl/add!
        o = ovec[n]
        # Set material parameter tensors and various object-related indices at Fw (= Uw or
        # Vw) locations.
        # The set parameter tensors are used in Maxwell's equations without modification if
        # subpixel smoothing is not needed at the locations lcmp.
        assign_param_obj!(gt, nw, param3d_cmp, obj3d_cmp, pind3d_cmp, oind3d_cmp, o, τlcmp)
    end

    return nothing
end

# Note that this function is used as a function barrier to achieve type stability of o.
function assign_param_obj!(gt::GridType,  # primal field (U) or dual field (V)
                           nw::Int,  # w = XX (1), YY (2), ZZ (3), grid node (4); not from 0 to 3, because 0 cannot be used as array index
                           param3d_cmp::AbsArr{CFloat,5},  # parameter array to set
                           obj3d_cmp::AbsArr{Object3,3},  # object array to set
                           pind3d_cmp::AbsArr{ParamInd,3},  # material parameter index array to set
                           oind3d_cmp::AbsArr{ObjInd,3},  # object index array to set
                           o::Object3,  # object whose parameter is taken
                           τlcmp::Tuple3{AbsVec{Float}})  # location of field components
    # Set the location indices of object boundaries.
    assert(all(issorted.(τlcmp)))
    bn, bp = bounds(o)  # (SVector{3}, SVector{3})
    subn = map((l,b) -> (n = findfirst(l.≥b); n==0 ? 1 : n), τlcmp, bn)  # SVec3Int
    subp = map((l,b) -> (n = findlast(l.≤b); n==0 ? length(l) : n), τlcmp, bp)  # SVec3Int
    I, J, K = map((nᵢ,nₑ) -> nᵢ:nₑ, subn, subp)  # SVector{3,UnitRange{Int}}

    # Assign param to param3d.
    gt′ = alter(gt)
    param = matparam(o,gt)
    pind′ = paramind(o,gt′)
    oind = objind(o)
    if o.shape isa Box{3,9} && (o.shape::Box{3,9}).p == @SMatrix(eye(3))  # o.shape is Cartesian box
        if nw == 4
            # Set the off-diagonal entries of the material parameter tensor.
            for nc = nXYZ, nr = next2(nc)  # column- and row-indices
                param3d_cmp[I, J, K, nr, nc] .= param[nr,nc]
            end
        else  # nw = 1, 2, 3
            # Set the diagonal entries of the material parameter tensor.
            param3d_cmp[I, J, K, nw, nw] .= param[nw,nw]
        end

        pind3d_cmp[I,J,K] .= pind′
        oind3d_cmp[I,J,K] .= oind
        obj3d_cmp[I,J,K] .= o
    else  # o.shape is not Cartesian box
        for k = K, j = J, i = I  # z-, y-, x-indices
            pt = t_ind(τlcmp, i, j, k)
            if pt ∈ o
                if nw == 4
                    # Set the off-diagonal entries of the material parameter tensor.
                    for nc = nXYZ, nr = next2(nc)  # column- and row-indices
                        param3d_cmp[i,j,k,nr,nc] = param[nr,nc]
                    end
                else  # nw = 1, 2, 3
                    # Set the diagonal entries of the material parameter tensor.
                    param3d_cmp[i,j,k,nw,nw] = param[nw,nw]
                end

                pind3d_cmp[i,j,k] = pind′
                oind3d_cmp[i,j,k] = oind
                obj3d_cmp[i,j,k] = o
            end  # if pt ∈ o
        end  # for k = ..., j = ..., i = ...
    end  # if o isa ...

    return nothing
end

function smooth_param!(param3d::Tuple2{AbsArr{CFloat,5}},  # parameter array to smooth
                       obj3d::Tuple24{AbsArr{<:Object3,3}},  # object array (does not change)
                       pind3d::Tuple24{AbsArr{ParamInd,3}},  # material parameter index array (does not chaneg)
                       oind3d::Tuple24{AbsArr{ObjInd,3}},  # object index array (does not change)
                       l::Tuple23{AbsVecReal},  # location of field components
                       l′::Tuple23{AbsVecReal},  # location of voxel corners without transformation by boundary conditions
                       σ::Tuple23{AbsVec{Bool}},  # false if on symmetry boundary
                       ∆τ′::Tuple23{AbsVecReal})  # amount of shift by BLOCH boundary conditions
    for ngt = nPD
        param3d_gt = param3d[ngt]
        gt = PD[ngt]
        for nw = 1:4  # w = XX, YY, ZZ, grid node
            # Set the grid types of the x-, y-, z-locations of Fw.
            gt_cmp = SVector(gt, gt, gt)
            gt_cmp = broadcast((k,w,g)->(k==w ? alter(g) : g), nXYZ, nw, gt_cmp)  # no change if nw = 4

            # Set the grid types of the voxel corners surroundnig Fw.
            gt_cmp′ = alter.(gt_cmp)

            # Choose vectors for Fw.
            lcmp = t_ind(l, gt_cmp)
            σcmp = t_ind(σ, gt_cmp)

            # Choose arrays and vectors for voxel corners.
            lcmp′ = t_ind(l′, gt_cmp′)
            ∆τcmp′ = t_ind(∆τ′, gt_cmp′)

            # Below, we don't use alter(ngt) because ngt components of obj3d, pind3d, oind3d
            # contain gt materials, even though they are defined at voxel corners.
            obj3d_cmp′ = obj3d[ngt][nw]
            pind3d_cmp′ = pind3d[ngt][nw]
            oind3d_cmp′ = oind3d[ngt][nw]

            # Set various arrays for the current component.
            smooth_param_cmp!(gt, nw, param3d_gt, obj3d_cmp′, pind3d_cmp′, oind3d_cmp′, lcmp, lcmp′, σcmp, ∆τcmp′)
        end
    end

    return nothing
end

# Below, XXX_cmp has size N, whereas XXX_cmp′ has size N+1 (and corresponds to voxel corners).
function smooth_param_cmp!(gt::GridType,  # primal field (U) or dual field (V)
                           nw::Int,  # w = XX (1), YY (2), ZZ (3), grid node (4)
                           param3d_gt::AbsArr{CFloat,5},  # parameter array to smooth
                           obj3d_cmp′::AbsArr{<:Object3,3},  # object array (does not change)
                           pind3d_cmp′::AbsArr{ParamInd,3},  # material parameter index array (does not chaneg)
                           oind3d_cmp′::AbsArr{ObjInd,3},  # object index array (does not change)
                           lcmp::Tuple3{AbsVecReal},  # location of field components
                           lcmp′::Tuple3{AbsVecReal},  # location of voxel corners without transformation by boundary conditions
                           σcmp::Tuple3{AbsVec{Bool}},  # false if on symmetry boundary
                           ∆τcmp′::Tuple3{AbsVecReal})  # amount of shift by BLOCH boundary conditions
    Nx, Ny, Nz = length.(lcmp)
    pind_vxl = Vector{Int}(8)  # material parameter indices inside voxel
    oind_vxl = Vector{Int}(8)  # object indices inside voxel
    ind_c = Vector{Int}(8)  # corner indices inside voxel

    for k = 1:Nz, j = 1:Ny, i = 1:Nx
        ijk_cmp = SVector(i,j,k)
        ijk_vxl = (ijk_cmp, ijk_cmp+1)  # Tuple2{SVec3Int}

        # Get objects at voxel corners.
        c = 0
        for kc = t_ind(ijk_vxl,nZ,nZ), jc = t_ind(ijk_vxl,nY,nY), ic = t_ind(ijk_vxl,nX,nX)
            c += 1
            pind_vxl[c] = pind3d_cmp′[ic,jc,kc]
            oind_vxl[c] = oind3d_cmp′[ic,jc,kc]
        end

        is_vxl_uniform = (pind_vxl[1]==pind_vxl[2]==pind_vxl[3]==pind_vxl[4]
                        ==pind_vxl[5]==pind_vxl[6]==pind_vxl[7]==pind_vxl[8])

        if !is_vxl_uniform  # smoothing only when nonuniform
            # Find Nparam_vxl (number of different material parameters inside the voxel).
            sort8!(ind_c, pind_vxl)  # pind_vxl[ind_c[n]] ≤ pind_vxl[ind_c[n+1]]
            Nparam_vxl, n_change = countdiff(ind_c, pind_vxl)  # n_change: last n where pind_vxl[ind_c[n]] changed

            # Depending on the value of Nparam_vxl, calculate param_cmp (averaged material
            # parameter for Fw).
            nout = @SVector zeros(3)
            rvol = 0.0
            if  Nparam_vxl == 2
                # Attempt to apply Kottke's subpixel smoothing algorithm.
                ind_c1, ind_c2 = ind_c[8], ind_c[1]::Int  # indices of corners with two material parameters
                sub_c1, sub_c2 = ind2sub((2,2,2), ind_c1), ind2sub((2,2,2), ind_c2)  # subscritpts of corners ind_c1 and ind_c2
                obj_c1, obj_c2 = obj3d_cmp′[t_ind(ijk_vxl, sub_c1).data...], obj3d_cmp′[t_ind(ijk_vxl, sub_c2).data...]

                with2objs = true
                oind_c1, oind_c2 = oind_vxl[ind_c1], oind_vxl[ind_c2]

                # Because the voxel has two different material parameters, it has either two
                # objects (when each parameter is composed of a single object), or more than
                # two objects (if either of the two parameters is composed of two or more
                # objects).

                # Check if param_c1 is composed of more than one object (such that the voxel
                # has more than two objects).
                for n = n_change:7  # n = 8 corresponds to oind_c1 and is omitted
                    (with2objs = oind_vxl[ind_c[n]]==oind_c1) || break
                end

                if with2objs  # single object for param_fg
                    # Check if param_c2 is composed of more than one object (such that the
                    # voxel has more than two objects).
                    for n = 2:n_change-1  # n = 1 corresponds oind_c2 and is omitted
                        (with2objs = oind_vxl[ind_c[n]]==oind_c2) || break
                    end
                end

                # Find which of obj_c1 and obj_c2 are the foreground and background objects.
                if oind_c1 > oind_c2
                    obj_fg, obj_bg = obj_c1, obj_c2
                    sub_fg = sub_c1
                else
                    assert(oind_c1≠oind_c2)
                    obj_fg, obj_bg = obj_c2, obj_c1
                    sub_fg = sub_c2
                end

                # Find
                # - param_fg (foreground material parameters),
                # - param_bg (background material parameters),
                # - nout (outward normal of the foreground object), and
                # - rvol (volume fraction of the foreground object inside the voxel).
                if !with2objs  # two material parameters but more than two objects in voxel
                    # In this case, the interface between two materials is not defined by
                    # the surface of a single object, so we have to estimate nout from the
                    # locations of the corners occupied by the two materials.
                    param_fg, param_bg, nout, rvol =
                        kottke_input_simple(ind_c, n_change, obj_fg, obj_bg)::Tuple{SMat3Complex,SMat3Complex,SVec3Float,Float}
                else  # Nobj_vxl == 2 (because Nobj_vxl ≠ 1 when Nparam_vxl == 2)
                    # When Nparam_vxl == Nobj_vxl == 2, different material parameters
                    # must correspond to different objects.
                    x₀ = t_ind(lcmp, ijk_cmp)  # SVec3Float: location of center of smoothing voxel
                    σvxl = t_ind(σcmp, ijk_cmp)
                    lvxl = t_ind(lcmp′, (ijk_cmp, ijk_cmp+1))
                    ∆fg = t_ind(∆τcmp′, ijk_cmp + SVector(sub_fg) - 1)  # SVec3Float; nonzero if corner ind_fg is outside periodic boundary

                    # See "Overall smoothing algorithm" above.
                    param_fg, param_bg, nout, rvol =
                        kottke_input_accurate(x₀, σvxl, lvxl, ∆fg, obj_fg, obj_bg, gt)::Tuple{SMat3Complex,SMat3Complex,SVec3Float,Float}
                end  # if Nobj_vxl

                if iszero(nout)  # includes case of Nparam_vxl ≥ 3
                    # Give up Kottke's subpixel smoothing and take simple averaging.
                    param_cmp = gt==PRIM ? amean_param(obj3d_cmp′, ijk_vxl, gt)::SMat3Complex :
                                           hmean_param(obj3d_cmp′, ijk_vxl, gt)::SMat3Complex
                else
                    # Perform Kottke's subpixel smoothing.
                    param_cmp = kottke_avg_param(param_fg, param_bg, nout, rvol)
                end
            end  # if Nparam_vxl == 2

            if nw == 4  # set off-diagonal entries of param using param_bg at grid nodes
                for nc = nXYZ, nr = next2(nc)  # column- and row-indices
                    param3d_gt[i, j, k, nr, nc] = param_cmp[nr,nc]
                end
            else  # w = x, y, z; set diagonal entries of param using param at Fw locations
                param3d_gt[i, j, k, nw, nw] = param_cmp[nw,nw]
            end
        end  # if Nparam_vxl ≠ 1
    end  # for kcmp = 1:Nz, jcmp = 1:Ny, icmp = 1:Nx

    return nothing
end

function amean_param(obj3d_cmp′::AbsArr{<:Object3,3}, ijk_vxl::Tuple2{SVec3Int}, gt::GridType)
    p = SMat3Complex(0,0,0, 0,0,0, 0,0,0)
    for kc = t_ind(ijk_vxl,nZ,nZ), jc = t_ind(ijk_vxl,nY,nY), ic = t_ind(ijk_vxl,nX,nX)
        o = obj3d_cmp′[ic,jc,kc]
        p += matparam(o.data.mat::EncodedMaterial, gt)
    end
    return p / length(ovec)
end

function hmean_param(obj3d_cmp′::AbsArr{<:Object3,3}, ijk_vxl::Tuple2{SVec3Int}, gt::GridType)
    p = SMat3Complex(0,0,0, 0,0,0, 0,0,0)
    for kc = t_ind(ijk_vxl,nZ,nZ), jc = t_ind(ijk_vxl,nY,nY), ic = t_ind(ijk_vxl,nX,nX)
        o = obj3d_cmp′[ic,jc,kc]
        p += inv(matparam(o.data.mat::EncodedMaterial, gt))
    end
    return inv(p / length(ovec))
end

function kottke_input_simple(ind_c::AbsVecInteger, n_change::Integer, obj_fg::Object3, obj_bg::Object3, gt::GridType)
    param_fg = matparam(obj_fg, gt)  # foreground material
    param_bg = matparam(obj_bg, gt)  # background material

    nout = @SVector zeros(3)  # nout for param_fg
    for n = n_change:8  # n = 8 corresponds ind_c[8] used for param_fg
        nout += NOUT_VXL[ind_c[n]]
    end
    rvol = (9-n_change) / 8

    return param_fg, param_bg, nout, rvol
end

function kottke_input_accurate(x₀::SVec3Float, σvxl::SVector{3,Bool}, lvxl::Tuple2{SVec3Float}, ∆fg::SVec3Float, obj_fg::Object3, obj_bg::Object3, gt::GridType)
    param_fg, param_bg = matparam(obj_fg, gt), matparam(obj_bg, gt)

    r₀, nout = surfpt_nearby(x₀ + ∆fg, obj_fg.shape)
    r₀ -= ∆fg
    nout = σvxl .* nout  # if voxel is across symmetry boundary plane, project nout to plane

    rvol = 0.0  # dummy value
    if !iszero(nout)
        rvol = volfrac(lvxl, nout, r₀)
    end

    return param_fg, param_bg, nout, rvol
end
