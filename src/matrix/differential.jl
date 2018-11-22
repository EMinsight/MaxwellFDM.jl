export create_∂, create_curl

# Assumption: we don't calculate derivatives for interpolated fields.  In other words, we
# calculate derivatives for the original fields defined on Yee's grid.  (There are
# exceptions; see below.)  This assumption simplifies implementation of the symmetry
# boundary condition in the derivative operators greatly.
#
# For example, suppose we use the primal grid planes as the E-field planes.  Then, for the
# symmetry boundary condition (which is the PEC boundary condition in this case), the fields
# defined on the domain boundaries are the tangential  E-fields and normal H-fields.  Both
# fields are zero on the symmetry boundaries.  Therefore, if we assume that we calculate
# derivatives only for those fields defined on Yee's grid, we do not have to construct
# separate derivative operators for the tangential E-field and normal H-field.  (Note that
# the derivative operator for the tangential E-field is used in curl for E, and the
# derivative operator for the normal H-field is used in divergence for H.)
#
# If we allow application of derivative operators to interpolated fields, there are cases
# where the boundary condition for those fields cannot be implemented in the derivative
# operators simply by zeroing the interpolated fileds.  For instance, in the aformentioned
# system, suppose we interpolate the normal E-field on the boundary.  Then, because the
# normal E-fields remain the same around the boundary, the interpolation on the boundary
# surface does not lead to zero.  Therefore, in the derivative operator for the normal
# E-field, we cannot implement the symmetry boundary condition simply by zeroing the fields.
# Supporting such interpolated fields would make the code for the derivative operators very
# convoluted, so we will assume that we don't deal with such cases.
#
# There are exceptions for this assumption.  As long as the interpolated fields on the
# symmetry boundaries are zero, our derivative operator code can be used to differentiate
# such interpolated fields.  For example, again in the aformentioned system, even though the
# normal D-field is not zero on the symmetry boundary, a portion of it could be zero in the
# following situation.  Suppose we are dealing with anisotropic materials.  Then, if we are
# dealing with the z-normal boundaries, Dz is the normal D-field but it is expressed as
# Dz = εzx Ex + εzy Ey + εzz Ez = Dzx + Dzy + Dzz.  Here, even though Dzz is nonzero because
# Ez is nonzero, Dzx and Dzy are zero because Ex and Ey are tangential to the PEC boundary.
# Therefore, we can use our derivative operator to differentiate Dzx and Dzy along the
# z-direction.  (I don't think there will be cases where we would need to differentiate Dzx
# and Dzy, but I mention this for consistency with the averaging operators in mean.jl,
# because we do need to take averages for Dzx and Dzy.)

## Discrete curl ##
create_curl(isfwd::AbsVecBool,  # isfwd[w] = true|false: create ∂w by forward|backward difference
            N::AbsVecInteger,  # size of grid
            ∆l::Tuple3{AbsVecNumber}=ones.((N...,)),  # ∆l[w]: distances between grid planes in x-direction
            isbloch::AbsVecBool=fill(true,length(N)),  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # Bloch phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    # I should not cast e⁻ⁱᵏᴸ into a complex vector, because then the entire curl matrix
    # becomes a complex matrix.  Sometimes I want to keep it real (e.g., when no PML and
    # Bloch phase factors are used).  In fact, this is the reason why I accept e⁻ⁱᵏᴸ instead
    # of constructing it from k and L as exp.(-im .* k .* L), which is always complex even
    # if k = 0.
    #
    # I should not cast ∆l to a vector of any specific type (e.g., Float, CFloat), either,
    # because sometimes I would want to even create an integral curl operator.
    (K = length(N); create_curl(SVector{K}(isfwd), SVector{K,Int}(N), ∆l, SVector{K}(isbloch), SVector{K}(e⁻ⁱᵏᴸ), reorder=reorder))


function create_curl(isfwd::SVec3Bool,  # isfwd[w] = true|false: create ∂w by forward|backward difference
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # ∆l[w]: distances between grid planes in x-direction
                     isbloch::SVec3Bool,  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVec3Number;  # Bloch phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    T = promote_type(eltype.(∆l)..., eltype(e⁻ⁱᵏᴸ))  # eltype(eltype(∆l)) can be Any if ∆l is inhomogeneous
    M = prod(N)

    Itot = VecInt(undef, 12M)
    Jtot = VecInt(undef, 12M)
    Vtot = Vector{T}(undef, 12M)

    indblk = 0  # index of matrix block
    for nv = nXYZ  # Cartesian compotent of output vector
        istr, ioff = reorder ? (3, nv-3) : (1, M*(nv-1))  # (row stride, row offset)
        parity = 1
        for nw = next2(nv)  # direction of differentiation
            nw′ = 6 - nv - nw  # Cantesian component of input vector; 6 = nX + nY + nZ
            jstr, joff = reorder ? (3, nw′-3) : (1, M*(nw′-1))  # (column stride, column offset)
            I, J, V = create_∂info(nw, isfwd[nw], N, ∆l[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw])

            @. I = istr * I + ioff
            @. J = jstr * J + joff
            V .*= parity

            # For some reason, using .= below is slower because it uses 1 allocatiotn.  On the
            # other hand, using = does not use allocation and therefore faster.
            indₛ, indₑ = indblk*2M + 1, (indblk+1)*2M
            Itot[indₛ:indₑ] = I
            Jtot[indₛ:indₑ] = J
            Vtot[indₛ:indₑ] = V
            indblk += 1

            parity = -1
        end
    end

    return dropzeros!(sparse(Itot, Jtot, Vtot, 3M, 3M))
end


## Difference operators ##
create_∂(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         isfwd::Bool,  # true|false for forward|backward difference
         N::AbsVecInteger,  # size of grid
         ∆w::Number=1.0,  # spatial discretization; vector of length N[nw]
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0) =  # Bloch phase factor
    (K = length(N); create_∂(nw, isfwd, N, fill(∆w, N[nw]), isbloch, e⁻ⁱᵏᴸ))  # fill: create vector of ∆w


create_∂(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         isfwd::Bool,  # true|false for forward|backward difference
         N::AbsVecInteger,  # size of grid
         ∆w::AbsVecNumber,  # spatial discretization; vector of length N[nw]
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0) =  # Bloch phase factor
    (K = length(N); M = prod(N); dropzeros!(sparse(create_∂info(nw, isfwd, SVector{K,Int}(N), ∆w, isbloch, e⁻ⁱᵏᴸ)..., M, M)))


# I need to figure out whether the ±1 entries of the backward difference operator is always
# the transpose of the forward difference operator for all boundary conditions. (∆w division
# factors are different, though.)  This was the case in the MATLAB code, but in the Julia
# code I changed the treatment of PDC, so let's make sure about this again.
#
# For the symmetry boundary for the E-field (then the E-fields are on the primal grid planes),
# assuming E that is correctly zeroed at the negative boundary is supplied to the forward
# difference operator, we can take the ±1 pattern of the forward difference operator for the
# Bloch boundary as the ±1 pattern of the forward difference operator for the symmetry
# boundary, because E = 0 is correctly used for the value at the positive boundary because
# the field value is copied from the negative boundary.
#
# Now, let's think about the backward difference operator for H.  Because the boundary is
# the symmetry boundary for E, the ghost H₀, which is before the negative boundary, must be
# the same as the non-ghost H₁.  Therefore, this leads to the first difference being H₁-H₀ =
# H₁-H₁ = 0, which means that the first row of the backward difference operator must be
# empty.  However, when the forward difference operator for the symmetry boundary is created
# the same as that for Bloch, its transpose does not have an empty first row!
#
# For the backward differce operater to be the transpose of the forward difference operator,
# I need to create the forward difference operator such that E₁ is zeroed.  This turns out
# to work.  See my notes on Sep/06/2017 in RN - MaxwellFDM.jl.

# Creates the w-directional difference matrix, with division by ∆w's.
function create_∂info(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
                      isfwd::Bool,  # true|false for forward|backward difference
                      N::SVector{K,Int},  # size of grid
                      ∆w::AbsVecNumber,  # spatial discretization; vector of length N[nw]
                      isbloch::Bool,  # boundary condition in w-direction
                      e⁻ⁱᵏᴸ::Number  # Bloch phase factor
                     ) where {K}
    M = prod(N)
    Nw = N[nw]
    ŵ = SVector(ntuple(identity,Val(K))) .== nw  # unit vector in w-direction; [0,true,0] for w == y
    ns = isfwd ? 1.0 : -1.0

    # Below, when constructing I, J's, V's, note that a tuple of array subscripts (i,j,k)
    # indicates a row index of the final matrix.  In other words, the entries with the same
    # subscripts correspond to the same *output* field.  So, when setting I[i,j,k], J₀[i,j,k],
    # Jₛ[i,j,k], V₀[i,j,k], Vₛ[i,j,k], we must know that they are for the entries in the
    # same row index (or the same output field).  Specifically,
    #
    # - I is the row index itself, so it is the identity map: I[i,j,k] = LinearIndices(N)[i,j,k].
    #
    # - J₀[i,j,k] is the column index of the diagonal entry in the row subscripted (i,j,k).
    # Because the row and column indices are the same for the diagonal entries, J₀ = I and
    # not created separately.
    #
    # - Jₛ[i,j,k] is the column index of the off-diagonal entry in the row subscripted
    # (i,j,k).
    #
    # - V₀[i,j,k] is the value of the diagonal entry in the row subscripted (i,j,k).
    #
    # - Vₛ[i,j,k] is the value fo the off-diagonal entry in the row subscripted (i,j,k).

    I = reshape(collect(1:M), N.data)

    Jₛ = reshape(collect(1:M), N.data)
    shifts = -ns * ŵ  # [0,-1,0] for w == y and ns = +1
    Jₛ = circshift(Jₛ, shifts.data)

    # Align ∆w in the w-direction.
    vec1 =  @SVector ones(Int,K)
    sizew = @. !ŵ * vec1 + ŵ * N  # [1,Ny,1] for w == y
    ∆W = reshape(∆w, sizew.data)

    # Construct the values of the diagonal and off-diagonal nonzero entries of the matrix
    # for the periodic (≠ Bloch) boundary condition.  (These values will be updated for
    # specific boundary conditions later.)
    #
    # To figure out the entries of V₀ and Vₛ, it is most convenient to figure out the entire
    # 3D arrays V₀ and Vₛ, instead of entry by entry.  In doing so, note that the entries
    # with the same subscripts correspond to the same *output* field (= same row index).
    T = promote_type(eltype(∆w), eltype(e⁻ⁱᵏᴸ))
    V₀ = -ns .* ones(T, N.data) ./ ∆W  # values of diagonal entries
    Vₛ = ns .* ones(T, N.data) ./ ∆W  # values of off-diagonal entries

    # Modify V₀ and Vₛ according to the boundary condition.  What we do is basically to take
    # the operator for the periodic boundary condition (constructed above) as a template and
    # then modify it for the actual boundary condition, either by applying the Bloch phase
    # factors for the Bloch boundary condition or by zeroing some diagonal and off-diagonal
    # entries for the symmetry boundary condition.  Because this modification is due to the
    # boundary condition, the modification needs to be made only at the array entries with
    # the first or last subscript of the Cartesion direction in which the differentiation is
    # taken.
    #
    # The outline of the procedure is as follows:
    #
    # 1. Detertmine the rows in which the modification should be made.  This is done by
    # figuring out the array subscripts (i,j,k) of the output derivative array entries to
    # which the input boundary fields contribute.
    #
    # 2. In each of these rows, determine the column in which the modification should be
    # made.  Because each row has only two nonzero entries (one diagoal and one off-diagonal),
    # this step is to determine whether the modification should be made at the diagonal
    # entry or the off-diagonal entry in the row currently concerned.  If the array
    # subscripts of the input fields contributing to the currently concerned output
    # derivatives are the same as the array subscripts of the output derivatives, the
    # modification is made in V₀ (diagonal).  If they are different, the modification is
    # made in Vₛ (off-diagonal).
    #
    # 3. In the determined rows, modify either V₀ or Vₛ depending on the above step.  This
    # is to modify V₀[i,j,k] or Vₛ[i,j,k], where (i,j,k) is the array subscripts specifying
    # the determined row.
    #
    #
    # For the final sparsity patterns of the operators, see my notes on September 6, 2017 in
    # RN - MaxwellFDM.jl.nb
    #
    # Below, Vₛ[Base.setindex(axes(Vₛ), iw, nw)...] mimics the implementation of slicedim
    # and basically means Vₛ[:,iw,:] for w = y.
    if isbloch
        # Application of the aformentioned procedure:
        #
        # 0. Goal
        # Our goal is to create the ghost fields on one boundary by multiplying the Bloch
        # phase factors to the nonghost fields on the opposite boundary.  Note that w is the
        # direction of differentiation.
        #
        # 1. Determination of the indices of the rows to modify
        # For isfwd = true, we are taking the forward difference between primal grid planes.
        # The ghost input fields are at the positive-end boundary, and they contritube to
        # the derivatives at the last dual grid points, so the rows to modify are
        # subscripted (i,j,k) with iw ≡ (i,j,k)[w] = Nw.
        #
        # For isfwd = false, we are taking the backward difference between dual grid planes.
        # The ghost input fields are before the negative-end boundary, and they contribute
        # to the derivatives at the first primal grid points, so the rows to modify are
        # subscripted (i,j,k) with iw ≡ (i,j,k)[w] = 1.
        #
        # 2. Determination of the column in each row to modify
        # For isfwd = true, the w-subscript of the ghost input field is greater by 1 than
        # that of the output derivative.
        #
        # For isfwd = false, the w-subscript of the ghost input field is less by 1 than that
        # of the output derivative.
        #
        # The nonghost input fields bing brought to these ghost fields (from the opposite
        # boundary) are different from the input fields subscripted the same as the output
        # derivatives.  Therefore, the matrix entries to modify by multiplying the Bloch
        # phase factors with are off-diagonal entries in Vₛ
        #
        # 3. Modification of the values in V
        # For isfwd = true, the nonghost input fields have a smaller w-subscript than the
        # ghost input fields that they are brought into.  Therefore, e⁻ⁱᵏᴸ must be
        # multiplied to the nonghost input fields to create the ghost input fields.
        #
        # For isfwd = false, the nonghost input fields have a larger w-subscript than the
        # ghost input fields that they are brought into.  Therefore, e⁺ⁱᵏᴸ must be
        # multiplied to the nonghost input fields to create the ghost input fields.
        iw = isfwd ? Nw : 1
        Vₛ[Base.setindex(axes(Vₛ), iw, nw)...] .*= e⁻ⁱᵏᴸ^ns
    else  # symmetry boundary
        # A. isfwd = true (forward difference)
        #
        # A-0. Goal
        # The forward difference is taken between primal grid planes.  Because the domain
        # boundaries are always primal grid planes, some input fields are on the domain
        # boundaries.  Those input fields are zero for the symmetry boundary condition and
        # do not contribute to the output derivatives.  Our goal is to make sure that the
        # forward derivative operator properly ignore those fields that should be zero.
        #
        # A-1. Determination of the indices of the the rows to modify
        # The negative-end boundary input fields contribute to the w-derivatives at the
        # first (nonghost) dual grid points in the w-direction.  The positive-end boundary
        # input fields contribute to the w-derivatives at the last dual grid points in the
        # w-direction.  Therefore, the array subscripts (i,j,k) we need to consider are the
        # ones with iw ≡ (i,j,k)[w] = 1 and Nw.
        #
        # A-2. Determination of the column in each row to modify
        # The negative-end boundary input fields have subscripts the same as the derivatives
        # they contribute to, so for iw = 1 we modify the diagonal entries stored in V₀.
        # The positive-end boundary input fields have w-subscript greater by 1 than the
        # derivatives they contribute to, so for iw = Nw we modify the off-diagonal entries
        # stored in Vₛ.
        #
        # A-3. Modification of the values in V
        # The negative-end boundary input fields must be zeroed, so for iw = 1 we set
        # V₀[i,j,k] = 0.
        # The positive-end boundary input fields must be also zeroed, so for iw = Nw we set
        # Vₛ[i,j,k] = 0.  Note that the positive-end boundary input fields are ghost fields
        # that are brought from the negative-end boundary (nonghost) input fields.
        #
        #
        # B. isfwd = false (backward difference)
        #
        # B-0. Goal
        # The backward difference is taken between dual grid planes.  The domain boundaries
        # are always primal grid planes.  Therefore, no input fields are on the domain
        # boundaries, and thus no input fields needs to be zeroed.
        #
        # However, the fields dual grid points are symmetric around the symmetry boundaries
        # (so that their interpolated values at the symmetry boundary are nonzero).
        # Therefore, the difference between the two is zero, leading to zero derivatives
        # evaluated at the symmetry boundaries.
        #
        # On the positive-end boundary, the derivative of the dual field is not evaluated,
        # so we don't have to do anything to realize the positivte-end symmetry boundary.
        # On the other hand, on the negative-end boundary the darivative of the dual field
        # is evaluated.  We have to make sure that the resulting derivatives are zero.
        #
        # B-1. Determination of the indices of the the rows to modify
        # The output derivatives at the negative-end boundary are subscripted (i,j,k) with
        # iw ≡ (i,j,k)[w] = 1.
        #
        # B-2. Determination of the column in each row to modify
        # To make the output derivatives at the negative-end boundary zero, both input
        # fields contributing to the derivatives must be zeroed.  Therefore, both V₀ and Vₛ
        # must be modified.
        #
        # B-3. Modification of the values in V
        # We simply zero the values, i.e., V₀[i,j,k] = Vₛ[i,j,k] = 0.
        #
        # The above procedure can be written compactly as follows without separating the
        # isfwd = true and false cases.


        # Zero the diagonal entries multiplied with the fields on the boundary.
        #
        # Note that the code below is independent of isfwd.  When the grid is uniorm, the
        # operators for the same boundary condition but for the opposite isfwd must be the
        # minus of the transpose of each other, so the diagonals of two operators must have
        # zeros at the same locations.  This means the locations of zeros on the diagonal
        # must be independent of isfwd.
        iw = 1
        V₀[Base.setindex(axes(V₀), iw, nw)...] .= 0

        # Zero the off-diagonal entries multiplied with the fields on the boundary.
        #
        # Note that the code below is indepent of the boundary condition.  When the grid is
        # uniform, this means that for the same ns, the operators for the symmetry boundary
        # conditions will have the same off-diagonal entries, except for the opposite signs.
        iw = isfwd ? Nw : 1
        Vₛ[Base.setindex(axes(Vₛ), iw, nw)...] .= 0
    end

    i = reshape(I, M)
    jₛ = reshape(Jₛ, M)
    v₀ = reshape(V₀, M)
    vₛ = reshape(Vₛ, M)

    Is = [i; i]  # row indices of [diagonal; off-diagonal]
    Js = [i; jₛ]  # column indices of [diagonal; off-diagonal]
    Vs = [v₀; vₛ]  # matrix entries of [diagonal, off-diagonal]

    return Is, Js, Vs
end
