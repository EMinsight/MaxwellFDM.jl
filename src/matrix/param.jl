export param3d2mat


param3d2mat(param3d::AbsArrComplex{5},
            gt::AbsVec{GridType},  # gt[w] = PRIM|DUAL for primal|dual grid planes as w-normal voxel faces
            N::AbsVecInteger,  # size of grid
            ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
            ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
            isbloch::AbsVecBool,  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # Bloch phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    (K = length(N); param3d2mat(param3d, SVector{K}(gt), SVector{K}(N), ∆l, ∆l′, SVector{K}(isbloch), SVector{K}(e⁻ⁱᵏᴸ), reorder=reorder))


function param3d2mat(param3d::AbsArrComplex{5},
                     gt::SVec3{GridType},  # gt[w] = PRIM|DUAL for primal|dual grid planes as w-normal voxel faces
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
                     ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
                     isbloch::SVec3Bool,  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVec3Number;  # Bloch phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    M = prod(N)

    # Following Oskooi et al.'s 2009 Optics Letters paper, off-diagonal entries of material
    # parameter tensors (e.g., ε) are evaluated at the corners of voxels whose edges are
    # field lines (e.g., E).  Therefore, if w-normal voxel faces are primal grid planes, the
    # input fields need to be averaged in the backward direction to be interpolated at voxel
    # corners.
    #
    # Once these interpolated input fields (e.g., Ex) are multiplied with off-diagonal
    # entries (e.g., εzx) of the material parameter tensor, we get the output fields (e.g.,
    # Dzx of Dz = Dzx + Dzy + Dzz = εzx Ex + εzy Ey + εzz Ez) in the direction normal (say
    # v-direction) to the input fields.  These output fields are still at the voxel corners.
    # Now, if the v-normal voxel faces are primal grid planes, these output fields need to
    # be averaged in the forward direction to be interpolated at voxel edges.
    #
    # In summary,
    #
    # - to obtain the fields to feed to the off-diagonal entries of material parameter
    # tensors, the w-component of the input fields need to be backward(forward)-averaged
    # along the w-direction if the w-normal voxel faces are primal (dual) grid planes, and
    #
    # - to distribute the resulting output fields in the v(≠w)-direction back to voxel edges,
    # the output fields need to be forward(backward)-averaged along the v-direction if the
    # v-normal voxel faces are primal (dual) grid planes.
    isfwd_in = gt.==Ref(DUAL)
    isfwd_out = gt.==Ref(PRIM)

    # For the output averaging, ∆l and ∆l′ are not supplied to create_mean in order to
    # create a simple arithmetic averaging operator.  This is because the area factor matrix
    # multiplied for symmetry to the left of the material parameter matrix multiplies the
    # same area factor to the two fields being averaged.  (See my notes on Jul/18/2018 in
    # MaxwellFDM in Agenda.)
    Mout = create_mean(isfwd_out, N, isbloch, e⁻ⁱᵏᴸ, reorder=reorder)


    kdiag = 0
    p3dmat = create_param3dmat(param3d, kdiag, N, reorder=reorder)  # diagonal components of ε tensor
    for kdiag = (1,-1)  # (superdiagonal, subdiagonal) components of ε tensor
        # For the input averaging, ∆l and ∆l′ are supplied to create min in order to create
        # a line integral averaging operator.  This is because the inverse of the length
        # factor matrix multiplied for symmetry to the right of the material parameter
        # matrix divides the two fields being averaged by different (= nonuniform) line
        # segments.  The ∆l factors multiplied inside create_minfo cancel the effect of this
        # multiplication with the nonuniform line segments.  (See my notes on Jul/18/2018 in
        # MaxwellFDM in Agenda.)
        Min = create_mean(isfwd_in, N, ∆l, ∆l′, isbloch, e⁻ⁱᵏᴸ, reorder=reorder)

        # Below, we use block-diagonal averaging matrices but either block-diagonal (kdiag
        # = 0), block-superdiagonal (kdiag = +1), or block-subdiagonal (kdiag = -1) material
        # parameter matrices.  See my bullet point entitled [Update (May/13/2018)] in
        # RN - Subpixel Smoothing.
        p3dmatₖ = create_param3dmat(param3d, kdiag, N, reorder=reorder)
        p3dmat += Mout * p3dmatₖ * Min
    end

    return p3dmat
end


function create_param3dmat(param3d::AbsArrComplex{5},
                           kdiag::Integer,  # 0|+1|-1 for diagonal|superdiagonal|subdiagonal of material parameter
                           N::SVec3Int;  # size of grid
                           reorder::Bool=true)  # true for more tightly banded matrix
    # Note that param3d's i, j, k indices run from 1 to N+1 rather than to N, so we should
    # not iterate those indices from 1 to end.
    M = prod(N)
    I = VecInt(undef, 3M)
    J = VecInt(undef, 3M)
    V = VecComplex(undef, 3M)
    n = 0
    for nv = nXYZ  # row index of tensor
        istr, ioff = reorder ? (3, nv-3) : (1, M*(nv-1))  # (row stride, row offset)
        nw = mod1(nv+kdiag, 3)
        jstr, joff = reorder ? (3, nw-3) : (1, M*(nw-1))  # (column stride, column offset)
        for k = 1:N[nZ], j = 1:N[nY], i = 1:N[nX]
            n += 1
            ind = LinearIndices(N.data)[i,j,k]  # linear index of Yee's cell

            I[n] = istr * ind + ioff
            J[n] = jstr * ind + joff
            V[n] = param3d[i,j,k,nv,nw]
        end
    end
    # for k = 1:N[nZ], j = 1:N[nY], i = 1:N[nX]
    #     ind = LinearIndices(N.data)[i,j,k]  # linear index of Yee's cell
    #     for nv = nXYZ  # row index of tensor
    #         n += 1
    #         istr, ioff = reorder ? (3, nv-3) : (1, M*(nv-1))  # (row stride, row offset)
    #         nw = mod1(nv+kdiag, 3)
    #         jstr, joff = reorder ? (3, nw-3) : (1, M*(nw-1))  # (column stride, column offset)
    #
    #         I[n] = istr * ind + ioff
    #         J[n] = jstr * ind + joff
    #         V[n] = param3d[i,j,k,nv,nw]
    #     end
    # end

    return sparse(I, J, V, 3M, 3M)
end
