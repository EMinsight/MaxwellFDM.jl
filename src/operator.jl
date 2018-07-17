export create_∂, create_curl, create_m, create_mean, create_param3dmat, param3d2mat

## Discrete curl ##
create_curl(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
            N::AbsVecInteger,  # size of grid
            ∆l::Tuple3{AbsVecNumber}=ones.((N...)),  # ∆l[w]: distances between grid planes in x-direction
            ebc::AbsVec{EBC}=fill(BLOCH,length(N)),  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # BLOCH phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    (K = length(N); create_curl(gt, SVector{K}(N), ∆l, SVector{K}(ebc), SVector{K}(e⁻ⁱᵏᴸ), reorder=reorder))


# I need to create create_curl_info! first.  Then, from there it is easy to eliminate some
# rows and columns from I, J, V.  I need to create a sparse matrix from such reduced I, J, V.
#
# Also in the future, change create_∂ to return only r, c, v vectors (instead of a sparse matrix)
# and create a sparse matrix at once.  This will create the curl matrix twice as fast.  I
# can even pre-permutate the collection of r's, c's, v's to create a permuted sparse matrix.
function create_curl(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # ∆l[w]: distances between grid planes in x-direction
                     ebc::SVector{3,EBC},  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVector{3,<:Number};  # BLOCH phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    ns = gt==PRIM ? 1 : -1
    T = promote_type(eltype.(∆l)..., eltype(e⁻ⁱᵏᴸ))  # eltype(eltype(∆l)) can be Any if ∆l is inhomogeneous
    M = prod(N)

    Itot = VecInt()
    Jtot = VecInt()
    Vtot = Vector{T}()

    for nv = nXYZ  # Cartesian compotent of output vector
        istr, ioff = reorder ? (3, nv-3) : (1, M*(nv-1))  # (row stride, row offset)
        parity = 1
        for nw = next2(nv)  # direction of differentiation
            nw′ = 6 - nv - nw  # Cantesian component of input vector; 6 = nX + nY + nZ
            jstr, joff = reorder ? (3, nw′-3) : (1, M*(nw′-1))  # (column stride, column offset)
            I, J, V = create_∂info(nw, ns, N, ∆l[nw], ebc[nw], e⁻ⁱᵏᴸ[nw])

            @. I = istr * I + ioff
            @. J = jstr * J + joff
            V .*= parity

            append!(Itot, I)
            append!(Jtot, J)
            append!(Vtot, V)

            parity = -1
        end
    end

    return sparse(Itot, Jtot, Vtot, 3M, 3M)
end


## Difference operators ##
create_∂(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ∆w::Number=1.0,  # spatial discretization; vector of length N[nw]
         ebc::EBC=BLOCH,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # BLOCH phase factor
        ) where {K} =
    create_∂(nw, ns, N, fill(∆w, N[nw]), ebc, e⁻ⁱᵏᴸ)  # fill: create vector of ∆w


create_∂(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ∆w::AbsVecNumber,  # spatial discretization; vector of length N[nw]
         ebc::EBC=BLOCH,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # BLOCH phase factor
        ) where {K} =
    (M = prod(N); sparse(create_∂info(nw, ns, N, ∆w, ebc, e⁻ⁱᵏᴸ)..., M, M))


# I need to figure out whether the ±1 entries of the backward difference operator is always
# the transpose of the forward difference operator for all boundary conditions. (∆w division
# factors are different, though.)  This was the case in the MATLAB code, but in the Julia
# code I changed the treatment of PDC, so let's make sure about this again.
#
# For PPC, assuming U with correctly zero at the negative boundary is supplied to the forward
# difference operator, the ±1 pattern of the difference operator must be the same as that of
# the BLOCH boundary, because 0 at the negative boundary is used for the value at the
# positive boundary.
#
# Now, let's think about the backward difference operator for V.  Because the boundary is
# PPC, the ghost V₀, which is before the negative boundary, must be the same as the
# non-ghost Vₛ.  Therefore, this leads to the first difference being Vₛ-V₀ = Vₛ-Vₛ = 0,
# which means that the first row of the backward difference operator must be empty.
# However, when the forward difference operator for PPC is created the same as that for
# BLOCH, its transpose does not have an empty first row!
#
# For the backward differce operater to be the transpose of the forward difference operator,
# I need to create the forward difference operator such that Uₛ is zeroed.  This turns out
# to work.  See the notes on Sep/06/2017 in RN - MaxwellFDM.jl.nb.

# Creates the w-directional difference matrix, with division by ∆w's.
function create_∂info(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
                      ns::Integer,  # 1|-1 for forward|backward difference
                      N::SVector{K,Int},  # size of grid
                      ∆w::AbsVecNumber,  # spatial discretization; vector of length N[nw]
                      ebc::EBC,  # boundary condition in w-direction
                      e⁻ⁱᵏᴸ::Number  # BLOCH phase factor
                     ) where {K}
    M = prod(N)
    Nw = N[nw]
    ŵ = SVector(ntuple(identity,Val{K})) .== nw  # unit vector in w-direction; [0,true,0] for w == y

    # Construct the row and column indices of nonzero entries of the matrix.
    I₀ = reshape(collect(1:M), N.data)  # row and column indices of diagonal entries
    Iₛ = reshape(collect(1:M), N.data)  # row indices of off-diagonal entries
    Jₛ = reshape(collect(1:M), N.data)  # column indices of off-diagonal entries
    shifts = -ns * ŵ  # [0,-1,0] for w == y and ns = +1
    Jₛ = circshift(Jₛ, shifts.data)

    # Align ∆w in the w-direction.
    vec1 =  @SVector ones(Int,K)
    sizew = @. !ŵ * vec1 + ŵ * N  # [1,Ny,1] for w == y
    ∆W = reshape(∆w, sizew.data)

    # Construct the values of the diagonal and off-diagonal nonzero entries of the matrix.
    T = promote_type(eltype(∆w), eltype(e⁻ⁱᵏᴸ))
    V₀ = -ns .* ones(T, N.data) ./ ∆W  # values of diagonal entries
    Vₛ = ns .* ones(T, N.data) ./ ∆W  # values of off-diagonal entries

    # Modify I, J, V according to the boundary condition; see my notes on September 6, 2017.
    if ebc == BLOCH
        iw = ns<0 ? 1 : Nw  # ghost points at negative end for ns < 0
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .*= e⁻ⁱᵏᴸ^ns  # mimic implementation of slicedim
    else  # ebc ≠ BLOCH
        # Set up the diagonal entries.  (This is independent of ns.)
        # Because the way to set up the diagonal entries doesn't depend on ns, the operators
        # for ns = +1 and –1 are the transpose of each other.
        iw = ebc==PPC ? 1 : Nw
        V₀[Base.setindex(indices(V₀), iw, nw)...] .= 0  # mimic implementation of slicedim

        # Set up the off-diagonal entries.  (This is indepent of boundary condition.)
        # The construction here guarantees the off-diagonal parts for ns = +1 and -1 are the
        # transpose of each other, regardless of boundary condition.
        iw = ns<0 ? 1 : Nw
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .= 0  # mimic implementation of slicedim
    end

    I = [I₀[:]; Iₛ[:]]  # row indices of [diagonal; off-diagonal]
    J = [I₀[:]; Jₛ[:]]  # column indices of [diagonal; off-diagonal]
    V = [V₀[:]; Vₛ[:]]  # matrix entries of [diagonal, off-diagonal]

    return I, J, V
end


# Creates the field-averaging operator for all three Cartegian components.
create_mean(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
            ns::Integer,  # 1|-1 for forward|backward averaging
            N::SVec3Int,  # size of grid
            ebc::SVector{3,EBC},  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::SVector{3,<:Number};  # BLOCH phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    create_mean(gt, ns, N, ones.(N.data), ones.(N.data), ebc, e⁻ⁱᵏᴸ, reorder=reorder)


function create_mean(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
                     ns::Integer,  # 1|-1 for forward|backward averaging
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
                     ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
                     ebc::SVector{3,EBC},  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVector{3,<:Number};  # BLOCH phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    T = promote_type(eltype.(∆l)..., eltype.(∆l′)..., eltype(e⁻ⁱᵏᴸ))  # eltype(eltype(∆l)) can be Any if ∆l is inhomogeneous
    M = prod(N)

    Itot = VecInt(6M)
    Jtot = VecInt(6M)
    Vtot = Vector{T}(6M)

    for nw = nXYZ  # Cartesian compotent of output vector
        indstr, indoff = reorder ? (3, nw-3) : (1, M*(nw-1))  # (row stride, row offset)

        I, J, V = create_minfo(gt, nw, ns, N, ∆l[nw], ∆l′[nw], ebc[nw], e⁻ⁱᵏᴸ[nw])

        @. I = indstr * I + indoff
        @. J = indstr * J + indoff

        indₛ, indₑ = 2(nw-1)*M+1, 2nw * M
        Itot[indₛ:indₑ] = I
        Jtot[indₛ:indₑ] = J
        Vtot[indₛ:indₑ] = V
    end

    return sparse(Itot, Jtot, Vtot, 3M, 3M)  # 3M×3M matrix with 2 entries per row (so 6M entries for I, J, V)
end


## Field-averaging operators ##
# This creates the averaging operator for a single Cartesian component.  For the operator
# for all three Cartesian components, use create_mean.
#
# Construction of these operators are similar to that of difference operators.  However,
# unlike the difference operators that are primarilly used for curl and hence differentiates
# the fields along the direction normal to the fields, the field-averaging operators average
# fields along the field direction.  As a result, backward (rather than forward) averaging
# is for primal fields.
create_m(gt::GridType,  # PRIM|DUAL for primal|dual field
         nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ebc::EBC=BLOCH,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # BLOCH phase factor
        ) where {K} =
    (M = prod(N); sparse(create_minfo(gt, nw, ns, N, ebc, e⁻ⁱᵏᴸ)..., M, M))


create_m(gt::GridType,  # PRIM|DUAL for primal|dual field
         nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ∆w::AbsVecNumber,  # line segments to multiply with; vector of length N[nw]
         ∆w′::AbsVecNumber,  # line segments to divide by; vector of length N[nw]
         ebc::EBC=BLOCH,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # BLOCH phase factor
        ) where {K} =
    (M = prod(N); sparse(create_minfo(gt, nw, ns, N, ∆w, ∆w′, ebc, e⁻ⁱᵏᴸ)..., M, M))


create_minfo(gt::GridType,  # PRIM|DUAL for primal|dual field
             nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
             ns::Integer,  # 1|-1 for forward|backward averaging
             N::SVector{K,Int},  # size of grid
             ebc::EBC,  # boundary condition in w-direction
             e⁻ⁱᵏᴸ::Number  # BLOCH phase factor
            ) where {K} =
    (∆w = ones(N[nw]); create_minfo(gt, nw, ns, N, ∆w, ∆w, ebc, e⁻ⁱᵏᴸ))


function create_minfo(gt::GridType,  # PRIM|DUAL for primal|dual field
                      nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
                      ns::Integer,  # 1|-1 for forward|backward averaging
                      N::SVector{K,Int},  # size of grid
                      ∆w::AbsVecNumber,  # line segments to multiply with; vector of length N[nw]
                      ∆w′::AbsVecNumber,  # line segments to divide by; vector of length N[nw]
                      ebc::EBC,  # boundary condition in w-direction
                      e⁻ⁱᵏᴸ::Number  # BLOCH phase factor
                     ) where {K}
    M = prod(N)
    Nw = N[nw]
    ŵ = SVector(ntuple(identity,Val{K})) .== nw  # [0,true,0] for w == y
    T = promote_type(eltype(∆w), eltype(∆w′), eltype(e⁻ⁱᵏᴸ))
    withcongbc = (gt==PRIM && ebc==PPC) || (gt==DUAL && ebc==PDC)  # bc type is congruent with field type

    # Arrange ∆w and ∆w′ in the w-direction by reshape.
    vec1 =  @SVector ones(Int,K)
    sizew = @. !ŵ * vec1 + ŵ * N  # [1,Ny,1] for w == y
    ∆W = reshape(∆w, sizew.data)
    ∆W′ = reshape(∆w′, sizew.data)

    # Construct the row indices and values of nonzero diagonal entries of the matrix.
    I₀ = reshape(collect(1:M), N.data)  # row and column indices of diagonal entries
    V₀ = fill(T(0.5), N.data) .* ∆W ./ ∆W′  # values of diagonal entries

    # Construct the row and column indices and values of nonzero off-diagonal entries of the
    # matrix.
    Iₛ = reshape(collect(1:M), N.data)  # row indices of off-diagonal entries
    Jₛ = reshape(collect(1:M), N.data)  # column indices of off-diagonal entries
    Vₛ = fill(T(0.5), N.data) .* ∆W  # values of off-diagonal entries (division later)

    shifts = -ns * ŵ  # [0,-1,0] for w == y and ns = +1
    Jₛ = circshift(Jₛ, shifts.data)
    Vₛ = circshift(Vₛ, shifts.data)
    Vₛ ./= ∆W′

    # Modify I, J, V according to the boundary condition; see my notes on September 6, 2017.
    if ebc == BLOCH
        iw = ns<0 ? 1 : Nw  # ghost points at negative end for ns < 0
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .*= e⁻ⁱᵏᴸ^ns  # mimic implementation of slicedim
    else  # ebc ≠ BLOCH
        # Set up the diagonal entries.  (This is independent of ns.)
        # Because the way to set up the diagonal entries doesn't depend on ns, the operators
        # for ns = +1 and –1 are the transpose of each other.

        # The part below is the same as create_∂info if withcongbc == false; otherwise, it
        # is different from create_∂info: 2 instead of 0 is used.  See my notes entitled
        # [Beginning of the part added on Sep/21/2017] in RN - Subpixel Smoothing.
        iw = ebc==PPC ? 1 : Nw
        val = withcongbc ? 2 : 0
        V₀[Base.setindex(indices(V₀), iw, nw)...] .*= val  # mimic implementation of slicedim

        # Set up the off-diagonal entries.  (This is indepent of boundary condition.)
        # The construction here guarantees the off-diagonal parts for ns = +1 and -1 are the
        # transpose of each other, regardless of boundary condition.
        iw = ns<0 ? 1 : Nw
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .= 0  # mimic implementation of slicedim
    end

    I = [I₀[:]; Iₛ[:]]  # row indices of [diagonal; off-diagonal]
    J = [I₀[:]; Jₛ[:]]  # column indices of [diagonal; off-diagonal]
    V = [V₀[:]; Vₛ[:]]  # matrix entries of [diagonal, off-diagonal]

    return I, J, V
end


param3d2mat(param3d::AbsArr{CFloat,5},
            gt::GridType,  # PRIM|DUAL for primal|dual field
            N::AbsVecInteger,  # size of grid
            ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
            ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
            ebc::AbsVec{EBC},  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # BLOCH phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    (K = length(N); param3d2mat(param3d, gt, SVector{K}(N), ∆l, ∆l′, SVector{K}(ebc), SVector{K}(e⁻ⁱᵏᴸ), reorder=reorder))


function param3d2mat(param3d::AbsArr{CFloat,5},
                     gt::GridType,  # PRIM|DUAL for primal|dual field
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
                     ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
                     ebc::SVector{3,EBC},  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVector{3,<:Number};  # BLOCH phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    M = prod(N)
    ns_in, ns_out = gt==PRIM ? (-1,1) : (1,-1)
    Mout = create_mean(gt, ns_out, N, ebc, e⁻ⁱᵏᴸ, reorder=reorder)

    kdiag = 0
    p3dmat = create_param3dmat(param3d, kdiag, N, reorder=reorder)  # diagonal components of ε tensor
    for kdiag = (1,-1)  # (superdiagonal, subdiagonal) components of ε tensor
        Min = create_mean(gt, ns_in, N, ∆l, ∆l′, ebc, e⁻ⁱᵏᴸ, reorder=reorder)
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
