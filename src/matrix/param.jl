export create_param3dmat, param3d2mat


param3d2mat(param3d::AbsArr{CFloat,5},
            gt::GridType,  # PRIM|DUAL for primal|dual field
            N::AbsVecInteger,  # size of grid
            ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
            ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
            isbloch::AbsVec{Bool},  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # Bloch phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    (K = length(N); param3d2mat(param3d, gt, SVector{K}(N), ∆l, ∆l′, SVector{K}(isbloch), SVector{K}(e⁻ⁱᵏᴸ), reorder=reorder))


function param3d2mat(param3d::AbsArr{CFloat,5},
                     gt::GridType,  # PRIM|DUAL for primal|dual field
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
                     ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
                     isbloch::SVec3Bool,  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVec3Number;  # Bloch phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    M = prod(N)

    # For primal fields, the input averaging is backward averaging (ns = -1) and the output
    # averaging is forward averaging (ns = +1).
    ns_in, ns_out = gt==PRIM ? (-1,1) : (1,-1)

    # For the output averaging, ∆l and ∆l′ are not supplied to create_mean in order to
    # create a simple arithmetic averaging operator.  This is because the area factor matrix
    # multiplied for symmetry to the left of the material parameter matrix multiplies the
    # same area factor to the two fields being averaged.  (See my notes on Jul/18/2018 in
    # MaxwellFDM in Agenda.)
    Mout = create_mean(gt, ns_out, N, isbloch, e⁻ⁱᵏᴸ, reorder=reorder)

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
        Min = create_mean(gt, ns_in, N, ∆l, ∆l′, isbloch, e⁻ⁱᵏᴸ, reorder=reorder)

        # Below, we use block-diagonal averaging matrices but either block-diagonal (kdiag
        # = 0), block-superdiagonal (kdiag = +1), or block-subdiagonal (kdiag = -1) material
        # parameter matrices.  See my bullet point entitled [Update (May/13/2018)] in
        # RN - Subpixel Smoothing.
        p3dmatₖ = create_param3dmat(param3d, kdiag, N, reorder=reorder)
        p3dmat += Mout * p3dmatₖ * Min
    end

    return p3dmat
end


function create_param3dmat(param3d::AbsArr{CFloat,5},
                           kdiag::Integer,  # 0|+1|-1 for diagonal|superdiagonal|subdiagonal of material parameter
                           N::SVec3Int;  # size of grid
                           reorder::Bool=true)  # true for more tightly banded matrix
    # Note that param3d's i, j, k indices run from 1 to N+1 rather than to N, so we should
    # not iterate those indices from 1 to end.
    M = prod(N)
    I = VecInt(3M)
    J = VecInt(3M)
    V = VecComplex(3M)
    n = 0
    for nv = nXYZ  # row index of tensor
        istr, ioff = reorder ? (3, nv-3) : (1, M*(nv-1))  # (row stride, row offset)
        nw = mod1(nv+kdiag, 3)
        jstr, joff = reorder ? (3, nw-3) : (1, M*(nw-1))  # (column stride, column offset)
        for k = 1:N[nZ], j = 1:N[nY], i = 1:N[nX]
            n += 1
            ind = sub2ind(N.data, i, j, k)  # linear index of Yee's cell

            I[n] = istr * ind + ioff
            J[n] = jstr * ind + joff
            V[n] = param3d[i,j,k,nv,nw]
        end
    end
    # for k = 1:N[nZ], j = 1:N[nY], i = 1:N[nX]
    #     ind = sub2ind(N.data, i, j, k)  # linear index of Yee's cell
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
