export create_m, create_mean

# Assumption: we don't calculate averages for interpolated fields.  In other words, we
# calculate averages for the original fields defined on Yee's grid.  (There are exceptions;
# see below.)  This assumption simplifies implementation of the symmetry boundary condition
# in the averaging operators greatly.
#
# For example, suppose we use the primal grid planes as the E-field planes.  Then, for the
# symmetry boundary condition (which is the PEC boundary condition in this case), the fields
# defined on the domain boundaries are the tangential  E-fields and normal H-fields.  Both
# fields are zero on the symmetry boundaries.  Therefore, if we assume that we calculate
# averages only for those fields defined on Yee's grid, we do not have to construct
# separate averaging operators for the tangential E-field and normal H-field.
#
# See the comments in differential.jl as well.

# The averaging scheme in SALT will be mostly used with Bloch, so the implementation of
# boundary condition is actually not that important.

# The symmetry of the averaged material parameters is guaranteed even though the forward
# and backward averaging operators are not symmetric when the symmetry boundary condition
# is used: one has element with the factor of 2 multiplied while the other with the factor
# of zero.  Still, the symmetry boundary enforces a special structure on the subpixel
# smoothed material parameter tensor on the symmetry boundary.  This tensor of a special
# structure has zero entries where the nonsymmetry between the forward and backward
# averaging matters and eventually makes the effect of this nonsymmetry does not propagate
# up to the symmetry of the final operator.  See my notes entitled [Beginning of the part
# added on Aug/14/2018] in RN - Subpixel Smoothing.nb.

# Creates the field-averaging operator for all three Cartegian components.
create_mean(isfwd::AbsVecBool,  # isfwd[w] = true|false for forward|backward averaging
            N::AbsVecInteger,  # size of grid
            isbloch::AbsVecBool=fill(true,length(N)),  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # Bloch phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    create_mean(isfwd, N, ones.(N.data), ones.(N.data), isbloch, e⁻ⁱᵏᴸ, reorder=reorder)

create_mean(isfwd::AbsVecBool,  # isfwd[w] = true|false for forward|backward averaging
            N::AbsVecInteger,  # size of grid
            ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
            ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
            isbloch::AbsVecBool=fill(true,length(N)),  # boundary conditions in x, y, z
            e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # Bloch phase factor in x, y, z
            reorder::Bool=true) =  # true for more tightly banded matrix
    (K = length(N); create_mean(SVector{K}(isfwd), SVector{K,Int}(N), ∆l, ∆l′, SVector{K}(isbloch), SVector{K}(e⁻ⁱᵏᴸ), reorder=reorder))

function create_mean(isfwd::SVec3Bool,  # isfwd[w] = true|false for forward|backward averaging
                     N::SVec3Int,  # size of grid
                     ∆l::Tuple3{AbsVecNumber},  # line segments to multiply with; vectors of length N
                     ∆l′::Tuple3{AbsVecNumber},  # line segments to divide by; vectors of length N
                     isbloch::SVec3Bool,  # boundary conditions in x, y, z
                     e⁻ⁱᵏᴸ::SVec3Number;  # Bloch phase factor in x, y, z
                     reorder::Bool=true)  # true for more tightly banded matrix
    T = promote_type(eltype.(∆l)..., eltype.(∆l′)..., eltype(e⁻ⁱᵏᴸ))  # eltype(eltype(∆l)) can be Any if ∆l is inhomogeneous
    M = prod(N)

    Itot = VecInt(undef, 6M)
    Jtot = VecInt(undef, 6M)
    Vtot = Vector{T}(undef, 6M)

    for nw = nXYZ  # Cartesian compotent of output vector
        indstr, indoff = reorder ? (3, nw-3) : (1, M*(nw-1))  # (row stride, row offset)

        I, J, V = create_minfo(nw, isfwd[nw], N, ∆l[nw], ∆l′[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw])

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
#
# This creates the averaging operator for a single Cartesian component.  For the operator
# for all three Cartesian components, use create_mean.
#
# Construction of these operators are similar to that of difference operators.  However,
# unlike the difference operators that are primarilly used for curl and hence differentiates
# the fields along the direction normal to the fields, the field-averaging operators average
# fields along the field direction.  As a result, backward (rather than forward) averaging
# is for primal fields.
create_m(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         isfwd::Bool,  # true|false for forward|backward averaging
         N::AbsVecInteger,  # size of grid
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0) =  # Bloch phase factor
    (K = length(N); ∆w = ones(N[nw]); create_m(nw, isfwd, SVector{K,Int}(N), ∆w, ∆w, isbloch, e⁻ⁱᵏᴸ))

create_m(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         isfwd::Bool,  # true|false for forward|backward averaging
         N::AbsVecInteger,  # size of grid
         ∆w::AbsVecNumber,  # line segments to multiply with; vector of length N[nw]
         ∆w′::AbsVecNumber,  # line segments to divide by; vector of length N[nw]
         isbloch::Bool=true,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0) =  # Bloch phase factor
    (K = length(N); M = prod(N); sparse(create_minfo(nw, isfwd, SVector{K,Int}(N), ∆w, ∆w′, isbloch, e⁻ⁱᵏᴸ)..., M, M))


function create_minfo(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
                      isfwd::Bool,  # true|false for forward|backward averaging
                      N::SVector{K,Int},  # size of grid
                      ∆w::AbsVecNumber,  # line segments to multiply with; vector of length N[nw]
                      ∆w′::AbsVecNumber,  # line segments to divide by; vector of length N[nw]
                      isbloch::Bool,  # boundary condition in w-direction
                      e⁻ⁱᵏᴸ::Number  # Bloch phase factor
                     ) where {K}
    M = prod(N)
    Nw = N[nw]
    ŵ = SVector(ntuple(identity,Val(K))) .== nw  # [0,true,0] for w == y
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

    # Align ∆w and ∆w′ in the w-direction.
    vec1 =  @SVector ones(Int,K)
    sizew = @. !ŵ * vec1 + ŵ * N  # [1,Ny,1] for w == y
    ∆W = reshape(∆w, sizew.data)
    ∆W′ = reshape(∆w′, sizew.data)

    # Construct the values of the diagonal and off-diagonal nonzero entries of the matrix
    # for the periodic (≠ Bloch) boundary condition.  (These values will be updated for
    # specific boundary conditions later.)
    #
    # To figure out the entries of V₀ and Vₛ, it is most convenient to figure out the entire
    # 3D arrays V₀ and Vₛ, instead of entry by entry.  In doing so, note that the entries
    # with the same subscripts correspond to the same *output* field (= same row index).
    T = promote_type(eltype(∆w), eltype(∆w′), eltype(e⁻ⁱᵏᴸ))
    V₀ = fill(T(0.5), N.data) .* ∆W ./ ∆W′  # values of diagonal entries

    Vₛ = fill(T(0.5), N.data) .* ∆W  # values of off-diagonal entries (division later)
    Vₛ = circshift(Vₛ, shifts.data)
    Vₛ ./= ∆W′

    # Modify V₀ and Vₛ according to the boundary condition.  What we do is basically to take
    # the operator for the periodic boundary condition (constructed above) as a template and
    # then modify it for the actual boundary condition, either by applying the Bloch phase
    # factors for the Bloch boundary condition or by multiplying 0 or 2 to some diagonal and
    # off-diagonal entries for the symmetry boundary condition.  Because this modification
    # is due to the boundary condition, the modification needs to be made only at the array
    # entries with the first or last subscript of the Cartesion direction in which averaging
    # is taken.
    #
    # The outline of the procedure is as follows:
    #
    # 1. Detertmine the rows in which the modification should be made.  This is done by
    # figuring out the array subscripts (i,j,k) of the output average array entry to which
    # the input boundary fields contribute.
    #
    # 2. In each of these rows, determine the column in which the modification should be
    # made.  Because each row has only two nonzero entries (one diagoal and one off-diagonal),
    # this step is to determine whether the modification should be made at the diagonal
    # entry or the off-diagonal entry in the row currently concerned.  If the array
    # subscripts of the input fields contributing to the currently concerned output averages
    # are the same as the array subscripts of the output averages, the modification is made
    # in V₀ (diagonal).  If they are different, the modification is made in Vₛ (off-diagonal).
    #
    # 3. In the determined rows, modify either V₀ or Vₛ depending on the above step.  This
    # is to modify V₀[i,j,k] or Vₛ[i,j,k], where (i,j,k) is the array subscripts specifying
    # the determined row.
    #
    # For the final sparsity patterns of the operators, see my notes entitled [Beginning of
    # the part added on Aug/14/2018] in RN - Subpixel Smoothing.nb.
    #
    # Below, Vₛ[Base.setindex(axes(Vₛ), iw, nw)...] mimics the implementation of slicedim
    # and basically means Vₛ[:,iw,:] for w = y.
    if isbloch
        # Application of the aformentioned procedure:
        #
        # 0. Goal
        # Our goal is to create the ghost fields on one boundary by multiplying the Bloch
        # phase factors to the nonghost fields on the opposite boundary.  Note that w is the
        # direction of averaging.
        #
        # 1. Determination of the indices of the rows to modify
        # For isfwd = true, we are taking forward averaging between primal grid planes.  The
        # ghost input fields are at the positive-end boundary, and they contritube to the
        # averages at the last dual grid points, so the rows to modify are subscripted
        # (i,j,k) with iw ≡ (i,j,k)[w] = Nw.
        #
        # For isfwd = false, we are taking backward averaging between dual grid planes.  The
        # ghost input fields are before the negative-end boundary, and they contribute to
        # the averages at the first primal grid points, so the rows to modify are
        # subscripted (i,j,k) with iw ≡ (i,j,k)[w] = 1.
        #
        # 2. Determination of the column in each row to modify
        # For isfwd = true, the w-subscript of the ghost input field is greater by 1 than
        # that of the output average.
        #
        # For isfwd = false, the w-subscript of the ghost input field is less by 1 than that
        # of the output average.
        #
        # The nonghost input fields being brought to these ghost fields (from the opposite
        # boundary) are different from the input fields subscripted the same as the output
        # averages.  Therefore, the matrix entries to modify by multiplying the Bloch phase
        # factors with are off-diagonal entries in Vₛ
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
    else  # symmetry bounndary
        # A. isfwd = true (forward averaging)
        #
        # A-0. Goal
        # Forward averaging is taken between primal grid planes.  Because the domain
        # boundaries are always primal grid planes, some input fields are on the domain
        # boundaries.  Those input fields are zero for the symmetry boundary condition and
        # do not contribute to the output averages.  Our goal is to make sure that the
        # forward averaging operator properly ignore those fields that should be zero.
        #
        # A-1. Determination of the indices of the the rows to modify
        # The negative-end boundary input fields contribute to the w-averages at the first
        # (nonghost) dual grid points in the w-direction.  The positive-end boundary input
        # fields contribute to the w-averages at the last dual grid points in the
        # w-direction.  Therefore, the array subscripts (i,j,k) we need to consider are the
        # ones with iw ≡ (i,j,k)[w] = 1 and Nw.
        #
        # A-2. Determination of the column in each row to modify
        # The negative-end boundary input fields have subscripts the same as the averages
        # they contribute to, so for iw = 1 we modify the diagonal entries stored in V₀.
        # The positive-end boundary input fields have w-subscript greater by 1 than the
        # averages they contribute to, so for iw = Nw we modify the off-diagonal entries
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
        # B. isfwd = false (backward averaging)
        #
        # B-0. Goal
        # Backward averaging is taken between dual grid planes.  The domain boundaries are
        # always primal grid planes.  Therefore, no input fields are on the domain
        # boundaries, and thus no input fields needs to be zeroed.
        #
        # However, the field at dual grid points are symmetric around the symmetry
        # boundaries (so that their interpolated values at the symmetry boundary are nonzero).
        # Therefore, the average between the two is double the contribution from one input
        # field, leading to multiplication of a factor of 2 at the symmetry boundaries.
        #
        # On the positive-end boundary, the average of the dual field is not evaluated,
        # so we don't have to do anything to realize the positivte-end symmetry boundary.
        # On the other hand, on the negative-end boundary the average of the dual field is
        # evaluated.  We have to make sure that the resulting average is double the
        # contribution from one input field there.
        #
        # B-1. Determination of the indices of the the rows to modify
        # The output averages at the negative-end boundary are subscripted (i,j,k) with
        # iw ≡ (i,j,k)[w] = 1.
        #
        # B-2. Determination of the column in each row to modify
        # The output averages at the negative-end boundary is taken between the first
        # (nonghost) dual grid point and the ghost point before that point.  We can realize
        # the symmetry boundary by doubling the effect of the input field at the first
        # dual grid point and zeroing the effect of the input field at the ghost point.
        #
        # The input field at the first dual grid point is subscripted the same as the output
        # average, so its effect is described in a diagonal entry in V₀.  The input field at
        # the ghost point is subscripted 1 less in the w-direction, so its effect is
        # described in an off-diagonal entry in Vₛ.
        #
        # B-3. Modification of the values in V
        # We double V₀[i,j,k] and zero Vₛ[i,j,k].
        #
        # The above procedure can be written compactly as follows without separating the
        # isfwd = true and false cases.


        # Replace some diagonal entries with 0 or 2.
        iw = 1
        val = isfwd ? 0 : 2
        V₀[Base.setindex(axes(V₀), iw, nw)...] .*= val

        # Replace some off-diagonal entries with 0.
        iw = isfwd ? Nw : 1
        Vₛ[Base.setindex(axes(Vₛ), iw, nw)...] .= 0

        # Note that the forward and backward averaging operators are not the transpose of
        # each other because they have different diagonals.  This may seem to make the final
        # material parameter matrix (with the input and output field averaging operators
        # multiplied) nonsymmetric.  However, zeros in the material parameter tensor
        # enforced by the symmetry boundary condition nullifies the possibility of
        # nonsymmetry.  See my notes entitled [Beginning of the part added on Aug/14/2018]
        # in RN - Subpixel Smoothing.nb.  See also kottke_input_accurate in smoothing.jl.
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
