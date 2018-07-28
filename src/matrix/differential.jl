export create_∂, create_curl

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
