# To-dos:
#
# - How do we handle TF/SF?  Need a capability to do subpixel smoothing only inside some
# region.
#
# - When an object touching a symmetry boundary is flipped into the space behind the
# boundary, the material axes must be flipped as well if the material is anisotropic.
# Currently, let's leave it unsupported.  In other words, anisotropic materials must touch
# the boundaries only when the boundary conditions are Bloch.

export smooth_param!

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

# Given a vector `v` of comparable objects and vector `ind` of indices that sort `v`, return
# the number `c` of different objects in `v` and the index `nₑ` of `ind` where the last new
# object occurs (i.e., v[ind[nₑ]] is the last new object).
#
# If v is uniform, nₑ = length(v) + 1 (i.e., it is out of bound).
#
# For efficiency, this function assumes `v` is a vector of objects assigned to a voxel, i.e.,
# length(v) = 8.
@inline function countdiff(ind::AbsVecInteger, v::AbsVec)
    c, nₑ = 1, 9
    @simd for n = 2:8
        @inbounds v[ind[n-1]] ≠ v[ind[n]] && (c += 1; nₑ = n)
    end
    return c, nₑ  # c: number of changes in v; nₑ: last n where v[ind[n]] changed
end

const NOUT_VXL =  # vectors from corners to center
[SVector(1.,1.,1.), SVector(-1.,1.,1.), SVector(1.,-1.,1.), SVector(-1.,-1.,1.),
SVector(1.,1.,-1.), SVector(-1.,1.,-1.), SVector(1.,-1.,-1.), SVector(-1.,-1.,-1.)]

# Overall smoothing algorithm
#
# Below, obj, pind, oind refers to an Object3 instance, parameter index (integer value) that
# distinguishes different material parameters, and oind that distinguishes different objects.
#
# - Assign obj, pind, oind to arrays object-by-object (see assignment.jl).
#     - For the locations of grid points to assign the objects to, use τlcmp (lcmp created considering BC).
# - Using pind and oind, determine voxels to perform subpixel smoothing.
# - Inside each of the voxels, figure out the foreground object with which subpixel smoothing is performed.
#     - Iterating voxel corners, find the object that is put latest.  That object is the foreground object.
#     - Complication occurs when the voxel corner assigned with the foreground object is outside the domain boundary.
#         - In that case, to calculate nout and rvol properly, we have to move the voxel center.  (This is effectively the same as moving the foreground object.)
#             - If the corner is outside the periodic boundary, translate the voxel center before surfpt_nearby (using ∆fg).
#             - If the corner is outside the symmetry boundary, zero the nout componnent normal to the boundary (using σvxl).
#             - ∆fg and σvxl can be obtained simply by looking at the indices and BC.

function smooth_param!(param3d::Tuple2{AbsArrComplex{5}},  # parameter array to smooth
                       obj3d::Tuple24{AbsArr{<:Object3,3}},  # object array (does not change)
                       pind3d::Tuple24{AbsArr{ParamInd,3}},  # material parameter index array (does not chaneg)
                       oind3d::Tuple24{AbsArr{ObjInd,3}},  # object index array (does not change)
                       l::Tuple23{AbsVecReal},  # location of field components
                       l′::Tuple23{AbsVecReal},  # location of voxel corners without transformation by boundary conditions
                       σ::Tuple23{AbsVecBool},  # false if on symmetry boundary
                       ∆τ′::Tuple23{AbsVecReal},  # amount of shift by Bloch boundary conditions
                       boundft::SVec3FT)  # boundary field type
    for nft = nEH
        param3d_ft = param3d[nft]
        ft = EH[nft]
        gt_cmp₀ = PD[2 .- (boundft .== ft)]  # grid type of voxel corners
        for nw = 1:4  # w = X̂, Ŷ, Ẑ, grid node
            # Set the grid types of the x-, y-, z-locations of Fw.
            gt_cmp = broadcast((k,w,g)->(k==w ? alter(g) : g), nXYZ, nw, gt_cmp₀)  # grid type of Fw; no change if nw = 4

            # Set the grid types of the voxel corners surroundnig Fw.
            gt_cmp′ = alter.(gt_cmp)

            # Choose vectors for Fw (which is at the centers of the voxel defined by the
            # voxel corners below).
            lcmp = t_ind(l, gt_cmp)
            σcmp = t_ind(σ, gt_cmp)

            # Choose arrays and vectors for voxel corners.
            lcmp′ = t_ind(l′, gt_cmp′)
            ∆τcmp′ = t_ind(∆τ′, gt_cmp′)

            # Below, we don't use alter(nft) because nft components of obj3d, pind3d, oind3d
            # contain ft materials.  (However, the locations where these materials were
            # evaluated are the complementary field locations, i.e., Hw locations for
            # electric materials if the present loop is smoothing the electric material
            # associated with Ew, because the Hw locations are the corners of the voxels
            # centered at the Ew locations.)
            obj3d_cmp′ = obj3d[nft][nw]
            pind3d_cmp′ = pind3d[nft][nw]
            oind3d_cmp′ = oind3d[nft][nw]

            # Set various arrays for the current component.
            smooth_param_cmp!(ft, nw, param3d_ft, obj3d_cmp′, pind3d_cmp′, oind3d_cmp′, lcmp, lcmp′, σcmp, ∆τcmp′)
        end
    end

    return nothing
end

# Below, XXX_cmp has size N, whereas XXX_cmp′ has size N+1 (and corresponds to voxel corners).
function smooth_param_cmp!(ft::FieldType,  # E- or H-field
                           nw::Int,  # w = X̂ (1), Ŷ (2), Ẑ (3), grid node (4)
                           param3d_ft::AbsArrComplex{5},  # parameter array to smooth
                           obj3d_cmp′::AbsArr{O,3},  # object array (does not change)
                           pind3d_cmp′::AbsArr{ParamInd,3},  # material parameter index array (does not chaneg)
                           oind3d_cmp′::AbsArr{ObjInd,3},  # object index array (does not change)
                           lcmp::Tuple3{AbsVecReal},  # location of field components
                           lcmp′::Tuple3{AbsVecReal},  # location of voxel corners without transformation by boundary conditions
                           σcmp::Tuple3{AbsVecBool},  # false if on symmetry boundary
                           ∆τcmp′::Tuple3{AbsVecReal}  # amount of shift by Bloch boundary conditions
                          ) where {O<:Object3}
    Nx, Ny, Nz = length.(lcmp)

    obj_vxl = Vector{O}(undef, 8)  # object inside voxel
    pind_vxl = VecInt(undef, 8)  # material parameter indices inside voxel
    oind_vxl = VecInt(undef, 8)  # object indices inside voxel
    ind_c = VecInt(undef, 8)  # corner indices inside voxel

    for k = 1:Nz, j = 1:Ny, i = 1:Nx
        ijk_cmp = SVector(i,j,k)
        ijk_vxl = (ijk_cmp, ijk_cmp.+1)  # Tuple2{SVec3Int}

        # Retrieve the elements assigned to voxel corners from 3D arrays.
        #
        # Unlike the rest of the code that is performed only where subpixel smoothing is
        # performed, the following for loop is performed at all grid points.  Therefore, the
        # number of assignments to perform in the following for loop can easily reach
        # several millions even for a simple 3D problem.  Simple assignment can
        # be surprisingly time-consuming when peformed such a large number of times.  The
        # use of @inbounds helps reducing assignment time significantly.
        c = 0
        for kc = t_ind(ijk_vxl,nZ,nZ), jc = t_ind(ijk_vxl,nY,nY), ic = t_ind(ijk_vxl,nX,nX)
            c += 1
            @inbounds obj_vxl[c] = obj3d_cmp′[ic,jc,kc]
            @inbounds pind_vxl[c] = pind3d_cmp′[ic,jc,kc]
            @inbounds oind_vxl[c] = oind3d_cmp′[ic,jc,kc]
        end

        # We could use Nparam_vxl calculated below instead of is_vxl_uniform: if Nparam_vxl
        # = 1, then is_vxl_uniform = true.  However, calculating Nparam_vxl requires calling
        # sort8! and countdiff, which are costly when called for all voxels.  Therefore,
        # by using is_vxl_uniform, we avoid calling those functions when unnecessary.
        is_vxl_uniform = (pind_vxl[1]==pind_vxl[2]==pind_vxl[3]==pind_vxl[4]
                        ==pind_vxl[5]==pind_vxl[6]==pind_vxl[7]==pind_vxl[8])

        if !is_vxl_uniform  # perform smoothing only when voxel is nonuniform
            # Find Nparam_vxl (number of different material parameters inside the voxel).
            sort8!(ind_c, pind_vxl)  # pind_vxl[ind_c[n]] ≤ pind_vxl[ind_c[n+1]]
            Nparam_vxl, n_change = countdiff(ind_c, pind_vxl)  # n_change: last n where pind_vxl[ind_c[n]] changed

            nout = @SVector zeros(3)
            rvol = 0.0

            # Calculate param_cmp (averaged material parameter for Fw) when the voxel is
            # composed of two material parameters (Nparam_vxl = 2).  This includes cases
            # where the voxel is composed of two material parameters but more than two
            # objects (i.e., more than one object in the voxel have the same material).
            #
            # If the voxel is composed of more than two material parameters (Nparam_vxl ≥ 3),
            # give up subpixel smoothing and leave the material parameters as originally
            # assigned.  This can happen, e.g., at the junction of three
            # materials.
            if  Nparam_vxl == 2
                # Attempt to apply Kottke's subpixel smoothing algorithm.
                with2objs = true  # even if Nparam_vxl = 2, there could be more than 2 objects, so will be updated
                ind_c1, ind_c2 = ind_c[8], ind_c[1]::Int  # indices of corners of two different material parameters
                obj_c1, obj_c2 = obj_vxl[ind_c1], obj_vxl[ind_c2]
                oind_c1, oind_c2 = oind_vxl[ind_c1], oind_vxl[ind_c2]

                # The corners ind_c[n_change:8] are occupied with one material parameter
                # because ind_c is sorted for pind_vxl.  However, those corners can still be
                # occupied with more than one object.  In that case, the voxel is composed
                # of more than two objects.
                #
                # Note that not two objects can logically mean one object or more than two
                # objects, but because Nparam_vxl ≥ 2, with2objs = false means more than two
                # objects inside the voxel.
                for n = 7:-1:n_change  # n = 8 corresponds to oind_c1 and is omitted
                    (with2objs = oind_vxl[ind_c[n]]==oind_c1) || break
                end

                # If corners ind_c[n_change:8] are composed of more than one objects, we
                # already know the voxel has more than two objects because the corners
                # ind_c[1:n_change-1] are guaranteed to have a different object than those
                # occupying the corners ind_c[n_change:8], so we don't have to test further
                # to see if the voxel has more than two objects.
                if with2objs  # single object for corners ind_c[n_change:8]
                    # The corners ind_c[1:n_change-1] are occupied with one material
                    # parameter because ind_c is sorted for pind_vxl, but those corners can
                    # still be occupied with more than one object.  In that case, the voxel
                    # is composed of more than two objects.
                    for n = 2:n_change-1  # n = 1 corresponds oind_c2 and is omitted
                        (with2objs = oind_vxl[ind_c[n]]==oind_c2) || break
                    end
                end

                # Find which of obj_c1 and obj_c2 are the foreground and background objects.
                #
                # When multiple objects have the same object index (because they are
                # essentially the same object across a periodic boundary), it doesn't matter
                # which object to choose as the foreground (or background) object, we
                # translate points properly across the domain (by ∆fg below) in order to
                # evaluate the surface normal direction on the correct object.  (At least
                # that is the intention, but this needs to be tested after assigning the
                # same object index is actually implemented.)
                if oind_c1 > oind_c2  # obj_c1 is foreground
                    obj_fg, obj_bg = obj_c1, obj_c2
                    ind_fg = ind_c1
                else  # obj_c2 is foreground
                    @assert oind_c1≠oind_c2
                    obj_fg, obj_bg = obj_c2, obj_c1
                    ind_fg = ind_c2
                end

                # Find
                # - param_fg (foreground material parameters),
                # - param_bg (background material parameters),
                # - nout (outward normal of the foreground object), and
                # - rvol (volume fraction of the foreground object inside the voxel).
                if !with2objs  # two material parameters but more than two objects in voxel
                    # In this case, the interface between two materials is not defined by
                    # the surface of a single object, so we estimate nout simply from the
                    # locations of the corners occupied by the two materials.
                    param_fg, param_bg, nout, rvol =
                        kottke_input_simple(ind_c, n_change, obj_fg, obj_bg, ft)::Tuple{SMat3Complex,SMat3Complex,SVec3Float,Float}
                else  # two objects
                    # When Nparam_vxl == Nobj_vxl == 2, different material parameters
                    # must correspond to different objects.
                    x₀ = t_ind(lcmp, ijk_cmp)  # SVec3Float: location of center of smoothing voxel
                    σvxl = t_ind(σcmp, ijk_cmp)
                    lvxl = (t_ind(lcmp′,ijk_cmp), t_ind(lcmp′,ijk_cmp.+1))

                    sub_fg = CartesianIndices((2,2,2))[ind_fg].I  # subscritpt of corner ind_fg
                    ∆fg = t_ind(∆τcmp′, ijk_cmp + SVector(sub_fg) .- 1)  # SVec3Float; nonzero if corner ind_fg is outside periodic boundary

                    # See "Overall smoothing algorithm" above.
                    param_fg, param_bg, nout, rvol =
                        kottke_input_accurate(x₀, σvxl, lvxl, ∆fg, obj_fg, obj_bg, ft)::Tuple{SMat3Complex,SMat3Complex,SVec3Float,Float}
                end  # if !with2objs

                if iszero(nout)  # includes case of Nparam_vxl ≥ 3
                    # Give up Kottke's subpixel smoothing and take simple averaging.
                    param_cmp = ft==EE ? amean_param(obj3d_cmp′, ijk_vxl, ft)::SMat3Complex :
                                           hmean_param(obj3d_cmp′, ijk_vxl, ft)::SMat3Complex
                else
                    # Perform Kottke's subpixel smoothing.
                    param_cmp = kottke_avg_param(param_fg, param_bg, nout, rvol)  # defined in material.jl
                end
            end  # if Nparam_vxl == 2

            if nw == 4  # set off-diagonal entries of param using param_bg at grid nodes
                for nc = nXYZ, nr = next2(nc)  # column- and row-indices
                    param3d_ft[i, j, k, nr, nc] = param_cmp[nr,nc]
                end
            else  # w = x, y, z; set diagonal entries of param using param at Fw locations
                param3d_ft[i, j, k, nw, nw] = param_cmp[nw,nw]
            end
        end  # if Nparam_vxl ≠ 1
    end  # for kcmp = 1:Nz, jcmp = 1:Ny, icmp = 1:Nx

    return nothing
end

function amean_param(obj3d_cmp′::AbsArr{<:Object3,3}, ijk_vxl::Tuple2{SVec3Int}, ft::FieldType)
    p = SMat3Complex(0,0,0, 0,0,0, 0,0,0)
    for kc = t_ind(ijk_vxl,nZ,nZ), jc = t_ind(ijk_vxl,nY,nY), ic = t_ind(ijk_vxl,nX,nX)
        o = obj3d_cmp′[ic,jc,kc]
        p += matparam(o, ft)
    end
    return p / 8
end

function hmean_param(obj3d_cmp′::AbsArr{<:Object3,3}, ijk_vxl::Tuple2{SVec3Int}, ft::FieldType)
    p = SMat3Complex(0,0,0, 0,0,0, 0,0,0)
    for kc = t_ind(ijk_vxl,nZ,nZ), jc = t_ind(ijk_vxl,nY,nY), ic = t_ind(ijk_vxl,nX,nX)
        o = obj3d_cmp′[ic,jc,kc]
        p += inv(matparam(o, ft))
    end
    return inv(p / 8)
end

function kottke_input_simple(ind_c::AbsVecInteger, n_change::Integer, obj_fg::Object3, obj_bg::Object3, ft::FieldType)
    param_fg = matparam(obj_fg, ft)  # foreground material
    param_bg = matparam(obj_bg, ft)  # background material

    nout = @SVector zeros(3)  # nout for param_fg
    for n = n_change:8  # n = 8 corresponds ind_c[8] used for param_fg
        nout += NOUT_VXL[ind_c[n]]
    end
    rvol = (9-n_change) / 8

    return param_fg, param_bg, nout, rvol
end

function kottke_input_accurate(x₀::SVec3Float, σvxl::SVec3Bool, lvxl::Tuple2{SVec3Float}, ∆fg::SVec3Float, obj_fg::Object3, obj_bg::Object3, ft::FieldType)
    param_fg, param_bg = matparam(obj_fg, ft), matparam(obj_bg, ft)

    r₀, nout = surfpt_nearby(x₀ + ∆fg, obj_fg.shape)
    r₀ -= ∆fg
    nout = σvxl .* nout  # if voxel is across symmetry boundary plane, project nout to plane

    rvol = 0.0  # dummy value
    if !iszero(nout)
        rvol = volfrac(lvxl, nout, r₀)
    end

    # The following block is to support the transformed anisotropic materials behind the
    # symmetry boundaries.  S * param * S is the material parameter behind the symmetry
    # boundary.  The code is trying to produce an averaged material parameter between the
    # original material and the symmetry material.
    #
    # The simple arithmetic averaging taken between the original and symmetry materials
    # below is not a very accurate treatment.  However, it produces the correct D-field
    # from the E-field on the symmetry boundary (where only the normal E-field exists).
    #
    # Having zeros at the right location in the material parameter tensor is critical for
    # achieving a symmetric matrix after field averaging!  See my notes entitled [Beginning
    # of the part added on Aug/14/2018] in RN - Subpixel Smoothing.nb.
    S = diagm(Val(0) => .!σvxl - σvxl)  # .!σvxl - σvxl = [1,-1,-1] for σvxl = [false,true,true] (x-normal symmetry boundary)
    param_fg = 0.5 * (param_fg + S * param_fg * S)
    param_bg = 0.5 * (param_bg + S * param_bg * S)

    return param_fg, param_bg, nout, rvol
end
