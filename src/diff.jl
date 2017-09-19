export create_∂, create_curl

create_∂(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ∆w::AbsVec{<:Number},  # spatial discretization; vector of length N[nw]
         ebc::EBC=BLOCH,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # BLOCH phase factor
        ) where {K} =
    (M = prod(N); sparse(create_∂info(nw, ns, N, ∆w, ebc, e⁻ⁱᵏᴸ)..., M, M))


create_∂(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
         ns::Integer,  # 1|-1 for forward|backward difference
         N::SVector{K,Int},  # size of grid
         ∆w::Number=1.0,  # spatial discretization; vector of length N[nw]
         ebc::EBC=BLOCH,  # boundary condition in w-direction
         e⁻ⁱᵏᴸ::Number=1.0  # BLOCH phase factor
        ) where {K} =
    create_∂(nw, ns, N, fill(∆w, N[nw]), ebc, e⁻ⁱᵏᴸ)  # fill: create vector of ∆w


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
# non-ghost V₁.  Therefore, this leads to the first difference being V₁-V₀ = V₁-V₁ = 0,
# which means that the first row of the backward difference operator must be empty.
# However, when the forward difference operator for PPC is created the same as that for
# BLOCH, its transpose does not have an empty first row!
#
# For the backward differce operater to be the transpose of the forward difference operator,
# I need to create the forward difference operator such that U₁ is zeroed.  This turns out
# to work.  See the notes on Sep/06/2017 in RN - MaxwellFD3D.jl.nb.

# Creates the w-directional difference matrix, with division by ∆w's.
function create_∂info(nw::Integer,  # 1|2|3 for x|y|z; 1|2 for horizontal|vertical
                      ns::Integer,  # 1|-1 for forward|backward difference
                      N::SVector{K,Int},  # size of grid
                      ∆w::AbsVec{<:Number},  # spatial discretization; vector of length N[nw]
                      ebc::EBC,  # boundary condition in w-direction
                      e⁻ⁱᵏᴸ::Number  # BLOCH phase factor
                     ) where {K}
    M = prod(N)
    Nw = N[nw]
    ŵ = SVector(ntuple(identity,Val{K})) .== nw  # [0,true,0] for w == YY

    # Construct the row and column indices of nonzero entries of the matrix.
    I₀ = reshape(collect(1:M), N.data)  # row and column indices of diagonal entries
    I₁ = reshape(collect(1:M), N.data)  # row indices of off-diagonal entries
    J₁ = reshape(collect(1:M), N.data)  # column indices of off-diagonal entries
    shifts = -ns * ŵ  # [0,-1,0] for w == YY and ns = +1
    J₁ = circshift(J₁, shifts.data)

    # Construct the values of the diagonal and off-diagonal nonzero entries of the matrix.
    vec1 =  @SVector ones(Int,K)
    wsize = @. !ŵ * vec1 + ŵ * N  # [1,Ny,1] for w == YY
    ∆W = reshape(∆w, wsize.data)  # align ∆w in the w-direction
    T = promote_type(eltype(∆w), eltype(e⁻ⁱᵏᴸ))
    V₀ = -ns .* ones(T, N.data) ./ ∆W  # values of diagonal entries
    V₁ = ns .* ones(T, N.data) ./ ∆W  # values of off-diagonal entries

    # Modify I, J, V according to the boundary condition; see my notes on September 6, 2017.
    if ebc == BLOCH
        if ns > 0
            # Ghost points are at the positive end.
            V₁[Base.setindex(indices(V₁), Nw, nw)...] .*= e⁻ⁱᵏᴸ  # mimic implementation of slicedim
        else  # ns < 0
            # Ghost points are at the negative end.
            V₁[Base.setindex(indices(V₁), 1, nw)...] .*= 1/e⁻ⁱᵏᴸ  # mimic implementation of slicedim
        end
    else  # ebc ≠ BLOCH
        # Diagonal entries
        if ebc == PPC
            I₀ = slicedim(I₀, nw, 2:Nw)
            V₀ = slicedim(V₀, nw, 2:Nw)
        else  # ebc == PDC
            I₀ = slicedim(I₀, nw, 1:Nw-1)
            V₀ = slicedim(V₀, nw, 1:Nw-1)
        end

        # Off-diagonal entries
        if ns > 0
            I₁ = slicedim(I₁, nw, 1:Nw-1)
            J₁ = slicedim(J₁, nw, 1:Nw-1)
            V₁ = slicedim(V₁, nw, 1:Nw-1)
        else  # ns < 0
            I₁ = slicedim(I₁, nw, 2:Nw)
            J₁ = slicedim(J₁, nw, 2:Nw)
            V₁ = slicedim(V₁, nw, 2:Nw)
        end
    end

    I = [I₀[:]; I₁[:]]  # row indices
    J = [I₀[:]; J₁[:]]  # column indices
    V = [V₀[:]; V₁[:]]  # matrix entries

    return I, J, V
end

create_curl(gt::GridType,
            N::AbsVec{<:Integer},
            ∆l::Tuple3{AbsVec{<:Number}},
            ebc::AbsVec{EBC},
            e⁻ⁱᵏᴸ::AbsVec{<:Number}=ones(length(N));
            reorder::Bool=true) =
    (K = length(N); create_curl(gt, SVector{K}(N), ∆l, SVector{K}(ebc), SVector{K}(e⁻ⁱᵏᴸ), reorder))


# I need to create create_curl_info! first.  Then, from there it is easy to eliminate some
# rows and columns from I, J, V.  I need to create a sparse matrix from such reduced I, J, V.
#
# Also in the future, change create_∂ to return only r, c, v vectors (instead of a sparse matrix)
# and create a sparse matrix at once.  This will create the curl matrix twice as fast.  I
# can even pre-permutate the collection of r's, c's, v's to create a permuted sparse matrix.
function create_curl(gt::GridType,
                     N::SVector{K,<:Integer},
                     ∆l::Tuple3{AbsVec{<:Number}},
                     ebc::SVector{K,EBC},
                     e⁻ⁱᵏᴸ::SVector{K,<:Number},
                     reorder::Bool
                    ) where {K}
    # Create the curl operator acting on the primal field U.

    # gt: PRIM for curl for primal fields, DUAL for curl for dual fields

    # Nx, Ny, Nz = N
    #
    # ∂yUz = create_∂(YY, (Nx+1,Ny+1,Nz), ∆l(Int(YY)))  # Nin = (Nx+1,Ny+1,Nz), Nout = (Nx+1,Ny,Nz)
    # ∂zUy = create_∂(ZZ, (Nx+1,Ny,Nz+1), ∆l(Int(ZZ)))  # Nin = (Nx+1,Ny,Nz+1), Nout = (Nx+1,Ny,Nz)
    # Zxx = spzeros(prod((Nx+1,Ny,Nz)), prod((Nx,Ny+1,Nz+1)))  # Nin = (Nx,Ny+1,Nz+1), Nout = (Nx+1,Ny,Nz)
    #
    # ∂zUx = create_∂(ZZ, (Nx,Ny+1,Nz+1), ∆l(Int(ZZ)))  # Nin = (Nx,Ny+1,Nz+1), Nout = (Nx,Ny+1,Nz)
    # ∂xUz = create_∂(XX, (Nx+1,Ny+1,Nz), ∆l(Int(XX)))  # Nin = (Nx+1,Ny+1,Nz), Nout = (Nx,Ny+1,Nz)
    # Zyy = spzeros(prod((Nx,Ny+1,Nz)), prod((Nx+1,Ny,Nz+1)))  # Nin = (Nx+1,Ny,Nz+1), Nout = (Nx,Ny+1,Nz)
    #
    # ∂xUy = create_∂(XX, (Nx+1,Ny,Nz+1), ∆l(Int(XX)))  # Nin = (Nx+1,Ny,Nz+1), Nout = (Nx,Ny,Nz+1)
    # ∂yUx = create_∂(YY, (Nx,Ny+1,Nz+1), ∆l(Int(YY)))  # Nin = (Nx,Ny+1,Nz+1), Nout = (Nx,Ny,Nz+1)
    # Zzz = spzeros(prod((Nx,Ny,Nz+1)), prod((Nx+1,Ny+1,Nz)))  # Nin = (Nx+1,Ny+1,Nz), Nout = (Nx,Ny,Nz+1)
    #
    # return [Zxx -∂zUy ∂yUz;
    #         ∂zUx Zyy -∂xUz;
    #         -∂yUx ∂xUy Zzz]

    ns = gt==PRIM ? 1 : -1
    T = promote_type(eltype.(∆l)..., eltype(e⁻ⁱᵏᴸ))  # eltype(eltype(∆l)) can be Any if ∆l is inhomogeneous
    M = prod(N)

    Itot = Vector{Int}()
    Jtot = Vector{Int}()
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

# function create_curlV(N::Tuple3{Int}, ∆l::Tuple3{AbsVec{T}}, ebc::Tuple32{EBC}, e⁻ⁱᵏᴸ::Tuple3{Number}=(1,1,1)) where {T<:Number}
#     # Create the curl operator acting on the dual field V.  Unlike the case for the primal
#     # field U, internal boundary conditions need to be supplied, because curlU can generate
#     # in-domain V from in-domain (including boundary) U, whereas curlV cannot generate
#     # in-domain U from in-domain V, but V needs to be expanded for ghost points.
#
#     !isproperebc(ebc, e⁻ⁱᵏᴸ) && throw(ArgumentError("ebc = $ebc and e⁻ⁱᵏᴸ = $e⁻ⁱᵏᴸ do not describe proper boundary conditions."))
#
#     Nin = (w -> N .+ (XYZ.==w)).(XYZ)
#     Nout = (w -> N .+ (XYZ.≠w)).(XYZ)
#     Min = prod.(Nin)
#     Mout = prod.(Nout)
#
#     curlV = spzeros(sum(Mout), sum(Min))
#     for r = XYZ
#         nr = Int(r)
#         outrange2 = 0
#         for k = 1:nr
#             outrange2 += Mout[k]
#         end
#         # outrange2 = sum(Mout[1:nr])  # Mout[1:nr] is type-unstable
#         outrange1 = outrange2 - Mout[nr] + 1
#         outrange = outrange1:outrange2
#         s = 1
#         for p = next3(r)
#             p == r && break
#             np = Int(p)
#             nq = 6 - (nr+np)
#             q = XYZ[nq]
#
#             isbloch = ebc[np][1]==BLOCH
#             factor = isbloch ? (e⁻ⁱᵏᴸ[np], 1/e⁻ⁱᵏᴸ[np]) : (-1).^(ebc[np].==PDC)
#             gapVq = create_ghostappender_dual(Nin[nq], p, isbloch, factor)
#             ∂pVq = create_∂(p, Nin[nq] .+ 2.*(XYZ.==p), ∆l[np])
#
#             inrange2 = 0
#             for k = 1:nq
#                 inrange2 += Min[k]
#             end
#             # inrange2 = sum(Min[1:nq])  # Min[1:nq] is type-unstable
#             inrange1 = inrange2 - Min[nq] + 1
#             inrange = inrange1:inrange2
#
#             curlV[outrange, inrange] = s * ∂pVq * gapVq
#             s = -1
#         end
#     end
#
#     return curlV
# end
