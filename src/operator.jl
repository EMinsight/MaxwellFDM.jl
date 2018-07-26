export create_∂, create_curl, create_m, create_mean, create_param3dmat, param3d2mat

## Discrete curl ##
create_curl(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
            N::AbsVecInteger,  # size of grid
            ∆l::Tuple3{AbsVecNumber}=ones.((N...)),  # ∆l[w]: distances between grid planes in x-direction
            isbloch::AbsVec{Bool}=fill(true,length(N)),  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # Bloch phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    # I should not cast e⁻ⁱᵏᴸ into a complex vector, because then the entire curl matrix
    # becomes a complex matrix.  Sometimes I want to keep it real (e.g., when no PML and
    # Bloch phase factors are used).
    #
    # I should not cast ∆l to a vector of any specific type (e.g., Float, CFloat), either,
    # because sometimes I would want to even create an integral curl operator.
    (K = length(N); create_curl(gt, SVector{K,Int}(N), ∆l, SVector{K,Bool}(isbloch), SVector{K}(e⁻ⁱᵏᴸ), reorder=reorder))


# I need to create create_curl_info! first.  Then, from there it is easy to eliminate some
# rows and columns from I, J, V.  I need to create a sparse matrix from such reduced I, J, V.
#
# Also in the future, change create_∂ to return only r, c, v vectors (instead of a sparse matrix)
# and create a sparse matrix at once.  This will create the curl matrix twice as fast.  I
# can even pre-permutate the collection of r's, c's, v's to create a permuted sparse matrix.
function create_curl(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # ∆l[w]: distances between grid planes in x-direction
                     isbloch::SVec3Bool,  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVec3Number;  # Bloch phase factor in x, y, z
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
            I, J, V = create_∂info(nw, ns, N, ∆l[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw])

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
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # Bloch phase factor
        ) where {K} =
    create_∂(nw, ns, N, fill(∆w, N[nw]), isbloch, e⁻ⁱᵏᴸ)  # fill: create vector of ∆w


create_∂(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ∆w::AbsVecNumber,  # spatial discretization; vector of length N[nw]
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # Bloch phase factor
        ) where {K} =
    (M = prod(N); sparse(create_∂info(nw, ns, N, ∆w, isbloch, e⁻ⁱᵏᴸ)..., M, M))


# I need to figure out whether the ±1 entries of the backward difference operator is always
# the transpose of the forward difference operator for all boundary conditions. (∆w division
# factors are different, though.)  This was the case in the MATLAB code, but in the Julia
# code I changed the treatment of PDC, so let's make sure about this again.
#
# For the symmetry boundary, assuming U with correctly zero at the negative boundary is
# supplied to the forward difference operator, the ±1 pattern of the difference operator
# must be the same as that of the Bloch boundary, because 0 at the negative boundary is used
# for the value at the positive boundary.
#
# Now, let's think about the backward difference operator for V.  Because the boundary is
# the symmetry boundary, the ghost V₀, which is before the negative boundary, must be the
# same as the non-ghost Vₛ.  Therefore, this leads to the first difference being Vₛ-V₀ =
# Vₛ-Vₛ = 0, which means that the first row of the backward difference operator must be
# empty.  However, when the forward difference operator for the symmetry boundary is created
# the same as that for Bloch, its transpose does not have an empty first row!
#
# For the backward differce operater to be the transpose of the forward difference operator,
# I need to create the forward difference operator such that Uₛ is zeroed.  This turns out
# to work.  See the notes on Sep/06/2017 in RN - MaxwellFDM.jl.

# Creates the w-directional difference matrix, with division by ∆w's.
function create_∂info(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
                      ns::Integer,  # 1|-1 for forward|backward difference
                      N::SVector{K,Int},  # size of grid
                      ∆w::AbsVecNumber,  # spatial discretization; vector of length N[nw]
                      isbloch::Bool,  # boundary condition in w-direction
                      e⁻ⁱᵏᴸ::Number  # Bloch phase factor
                     ) where {K}
    M = prod(N)
    Nw = N[nw]
    ŵ = SVector(ntuple(identity,Val{K})) .== nw  # unit vector in w-direction; [0,true,0] for w == y

    # Construct the row and column indices of nonzero entries of the matrix.
    I₀ = reshape(collect(1:M), N.data)  # row and column indices of diagonal entries; row indices of off-diagonal entries
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

    # Modify I, J, V according to the boundary condition.  What we do is basically to take
    # the operator for Bloch as a template and then modify it for the symmetry boundary by
    # zeroing some diagonal and off-diagonal entries.
    #
    # Here are the details of what we do for ns = +1 (forward differentiation).  Note that
    # the operators for ns = +1 are applied to primal fields parallel to the boundaries.
    #
    # our goal is to make sure the primal fields on the boundaries do not contribute to the
    # derivatives at the first (nonghost) and last dual grid points.  For the derivatives at
    # the first dual grid points, the goal is achieved by zeroing the Bloch operator's
    # (diagonal) entries multiplied with the primal fields defined on the negative end.
    # For the derivatives at the last dual grid points, the goal is achieved by zeroing the
    # Bloch operator's off-diagonal entries multiplied with the same fields as above,
    # because these fields are periodically wrapped around and used as the (ghost) primal
    # fields defined on the positive end.
    #
    # For the final sparsity patterns of the operators, see my notes on September 6, 2017 in
    # RN - MaxwellFDM.jl.
    #
    # Below, Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] mimics the implementation of slicedim
    # and basically means Vₛ[:,iw,:] for w = y.
    if isbloch
        # - For ns = +1, multiply the negative-end field with e⁻ⁱᵏᴸ to bring it to the ghost
        # point at the positive end.
        # - For ns = -1, divide the positive-end field by e⁻ⁱᵏᴸ to bring it to the ghost
        # point at the negative end.
        iw = ns<0 ? 1 : Nw  # ghost points at negative end for ns < 0
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .*= e⁻ⁱᵏᴸ^ns
    else  # symmetry boundary
        # Zero the diagonal entries multiplied with the fields on the boundary.
        #
        # Note that the code below is independent of ns.  When the grid is uniorm, the
        # operators for the same boundary condition but for the opposite ns' must be the
        # minus of the transpose of each other, so the diagonals of two operators must have
        # zeros at the same locations.  This means the locations of zeros on the diagonal
        # must be independent of ns.
        iw = 1
        V₀[Base.setindex(indices(V₀), iw, nw)...] .= 0

        # Zero the off-diagonal entries multiplied with the fields on the boundary.
        #
        # Note that the code below is indepent of the boundary condition.  When the grid is
        # uniform, this means that for the same ns, the operators for the symmetry boundary
        # conditions will have the same off-diagonal entries, except for the opposite signs.
        iw = ns<0 ? 1 : Nw
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .= 0
    end

    i₀ = reshape(I₀, M)
    jₛ = reshape(Jₛ, M)
    v₀ = reshape(V₀, M)
    vₛ = reshape(Vₛ, M)

    I = [i₀; i₀]  # row indices of [diagonal; off-diagonal]
    J = [i₀; jₛ]  # column indices of [diagonal; off-diagonal]
    V = [v₀; vₛ]  # matrix entries of [diagonal, off-diagonal]

    return I, J, V
end


# Creates the field-averaging operator for all three Cartegian components.
create_mean(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
            ns::Integer,  # 1|-1 for forward|backward averaging
            N::SVec3Int,  # size of grid
            isbloch::SVec3Bool,  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::SVec3Number;  # Bloch phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    create_mean(gt, ns, N, ones.(N.data), ones.(N.data), isbloch, e⁻ⁱᵏᴸ, reorder=reorder)


function create_mean(gt::GridType,  # PRIM|DUAL for curl on primal|dual grid
                     ns::Integer,  # 1|-1 for forward|backward averaging
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
                     ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
                     isbloch::SVec3Bool,  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVec3Number;  # Bloch phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    T = promote_type(eltype.(∆l)..., eltype.(∆l′)..., eltype(e⁻ⁱᵏᴸ))  # eltype(eltype(∆l)) can be Any if ∆l is inhomogeneous
    M = prod(N)

    Itot = VecInt(6M)
    Jtot = VecInt(6M)
    Vtot = Vector{T}(6M)

    for nw = nXYZ  # Cartesian compotent of output vector
        indstr, indoff = reorder ? (3, nw-3) : (1, M*(nw-1))  # (row stride, row offset)

        I, J, V = create_minfo(gt, nw, ns, N, ∆l[nw], ∆l′[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw])

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
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # Bloch phase factor
        ) where {K} =
    (M = prod(N); sparse(create_minfo(gt, nw, ns, N, isbloch, e⁻ⁱᵏᴸ)..., M, M))


create_m(gt::GridType,  # PRIM|DUAL for primal|dual field
         nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ∆w::AbsVecNumber,  # line segments to multiply with; vector of length N[nw]
         ∆w′::AbsVecNumber,  # line segments to divide by; vector of length N[nw]
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # Bloch phase factor
        ) where {K} =
    (M = prod(N); sparse(create_minfo(gt, nw, ns, N, ∆w, ∆w′, isbloch, e⁻ⁱᵏᴸ)..., M, M))


create_minfo(gt::GridType,  # PRIM|DUAL for primal|dual field
             nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
             ns::Integer,  # 1|-1 for forward|backward averaging
             N::SVector{K,Int},  # size of grid
             isbloch::Bool,  # boundary condition in w-direction
             e⁻ⁱᵏᴸ::Number  # Bloch phase factor
            ) where {K} =
    (∆w = ones(N[nw]); create_minfo(gt, nw, ns, N, ∆w, ∆w, isbloch, e⁻ⁱᵏᴸ))


function create_minfo(gt::GridType,  # PRIM|DUAL for primal|dual field
                      nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
                      ns::Integer,  # 1|-1 for forward|backward averaging
                      N::SVector{K,Int},  # size of grid
                      ∆w::AbsVecNumber,  # line segments to multiply with; vector of length N[nw]
                      ∆w′::AbsVecNumber,  # line segments to divide by; vector of length N[nw]
                      isbloch::Bool,  # boundary condition in w-direction
                      e⁻ⁱᵏᴸ::Number  # Bloch phase factor
                     ) where {K}
    M = prod(N)
    Nw = N[nw]
    ŵ = SVector(ntuple(identity,Val{K})) .== nw  # [0,true,0] for w == y
    T = promote_type(eltype(∆w), eltype(∆w′), eltype(e⁻ⁱᵏᴸ))
    withcongbc = (gt==PRIM && !isbloch)  # bc type is congruent with field type

    # Arrange ∆w and ∆w′ in the w-direction by reshape.
    vec1 =  @SVector ones(Int,K)
    sizew = @. !ŵ * vec1 + ŵ * N  # [1,Ny,1] for w == y
    ∆W = reshape(∆w, sizew.data)
    ∆W′ = reshape(∆w′, sizew.data)

    # Construct the row indices and values of nonzero diagonal entries of the matrix.
    I₀ = reshape(collect(1:M), N.data)  # row and column indices of diagonal entries; row indices of off-diagonal entries
    V₀ = fill(T(0.5), N.data) .* ∆W ./ ∆W′  # values of diagonal entries

    # Construct the row and column indices and values of nonzero off-diagonal entries of the
    # matrix.
    Jₛ = reshape(collect(1:M), N.data)  # column indices of off-diagonal entries
    Vₛ = fill(T(0.5), N.data) .* ∆W  # values of off-diagonal entries (division later)

    shifts = -ns * ŵ  # [0,-1,0] for w == y and ns = +1
    Jₛ = circshift(Jₛ, shifts.data)
    Vₛ = circshift(Vₛ, shifts.data)
    Vₛ ./= ∆W′

    # Modify I, J, V according to the boundary condition.  What we do is basically to take
    # the operator for Bloch as a template and then modify it for the symmetry boundary by
    # replacing ±1's with 0's and 2's.
    #
    # Note that the operators for ns = -1 (backward averaging) and ns = +1 (forward
    # averaging) can be created for both the primal and dual fields.  (That is why
    # create_minfo takes an additional argument gt::GridType that specifies the type of the
    # field compared to create_∂info.)  The difference between the operators for ns = ±1 is
    # that one is for averaging the regular primal or dual fields, whereas the other is for
    # averaging the primal or dual fields already averaged at their corresponding voxel
    # corners.  Which of ns = +1 and -1 correspond to the regular and already-averaged
    # fields depends on whether the fields are primal or dual.  (This is easy to figure out.)
    #
    # Here are the details of what we do for ns = -1 (backward averaging) for primal fields.
    #
    # Our goal is to make sure symmetry around the symmetry boundary is correctly accounted
    # for during averaging.
    # For the averages at the primal grid point on the negative-end boundary, the goal is
    # achieved by substituting 2's for the Bloch operator's (diagonal) 1's that are
    # multiplied with the first (nonghost) primal fields (defined on the dual grid points),
    # and 0's for the Bloch operator's (off-diagonal) 1's that are multiplied with the ghost
    # primal fields, which are copied from the last (nonghost) primal fields.
    # For the averages at the primal grid point on the positive-end boundary, we don't have
    # to do anything because these primal grid points are ghost points and therefore the
    # averages are neither evaluated nor requested there.  (For example, the averaging of Ex
    # in the x-direction is performed to find the contribution of Ex to Dz and Dy, but Dz
    # and Dy on the ghost primal grid points on the positive-end boundary are ghost fields,
    # so the contribution of Ex to them are not requested during the equation construction.)
    #
    # For ns = +1 (forward averaging) for primal fields, it turns out that we only need to
    # replace entries with 0's, not with 2's.  This causes an asymmetry issue between the
    # ns = ±1 operators.  (For example, if the y-boundaries are the symmetry boundaries, the
    # backward averaging operators with 2's are multiplied to the right of ξxy to average
    # the input Uy fields, but the forward averaging operators without 2's are multiplied to
    # the left of ξ_yx to average the output [ξU]y fields.  Even if ξxy = ξyx, the resulting
    # material parameter matrix is not symmetric.)  However, it turns out that we can still
    # use 2's in the forward averaging operators to recover symmetry, because the fields to
    # be multiplied with these "wrong" 2's in the forward avegaring operator are guaranteed
    # to be zero. (In the above example, the output [ξU]y fields are created by ξyx Ux and
    # ξyz Uz, but Ux and Uz are zeros on the symmetry y-boundaries, so the [ξU]y to average
    # is already zero.)
    #
    # For more details and the final sparsity patterns of the operators, see my notes
    # entitled [Beginning of the part added on Sep/21/2017] in RN - Subpixel Smoothing.  See
    # also the notes on September 6, 2017 in RN - MaxwellFDM.jl for the matrices created by
    # create_∂info, because overall similar principles are used.
    #
    # Below, Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] mimics the implementation of slicedim
    # and basically means Vₛ[:,iw,:] for w = y.
    if isbloch
        # - For ns = +1, multiply the negative-end field with e⁻ⁱᵏᴸ to bring it to the ghost
        # point at the positive end.
        # - For ns = -1, divide the positive-end field by e⁻ⁱᵏᴸ to bring it to the ghost
        # point at the negative end.
        iw = ns<0 ? 1 : Nw  # ghost points at negative end for ns < 0
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .*= e⁻ⁱᵏᴸ^ns
    else  # symmetry bounndary
        # Replace some diagonal entries with 0 or 2.
        #
        # Note that the code below is independent of ns.  When the grid is uniorm, the
        # operators for the same boundary condition but for the opposite ns' must be the
        # minus of the transpose of each other, so the diagonals of two operators must have
        # zeros at the same locations.  This means the locations of zeros on the diagonal
        # must be independent of ns.
        #
        # The part below is the same as create_∂info if withcongbc == false; otherwise, it
        # is different from create_∂info: 2 instead of 0 is used.  See my notes entitled
        # [Beginning of the part added on Sep/21/2017] in RN - Subpixel Smoothing.
        iw = isbloch ? Nw : 1
        val = withcongbc ? 2 : 0
        V₀[Base.setindex(indices(V₀), iw, nw)...] .*= val

        # Replace some off-diagonal entries with 0.  Note that the code below is indepent of
        # the boundary condition.
        iw = ns<0 ? 1 : Nw
        Vₛ[Base.setindex(indices(Vₛ), iw, nw)...] .= 0
    end

    i₀ = reshape(I₀, M)
    jₛ = reshape(Jₛ, M)
    v₀ = reshape(V₀, M)
    vₛ = reshape(Vₛ, M)

    I = [i₀; i₀]  # row indices of [diagonal; off-diagonal]
    J = [i₀; jₛ]  # column indices of [diagonal; off-diagonal]
    V = [v₀; vₛ]  # matrix entries of [diagonal, off-diagonal]

    return I, J, V
end


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
