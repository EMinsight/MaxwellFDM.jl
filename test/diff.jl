@testset "differential" begin

@testset "create_∂" begin
    Nmax = 10
    N = SVector(rand(1:Nmax, 3)...)
    # N = SVector(3,3,3)
    M = prod(N)
    for nw = nXYZ
    # for nw = (1,)
        Nw = N[nw]
        sub′ = Vector{Int}(3)

        for ns = (-1,1)
        # for ns = (1,)
            ∂ws = spzeros(M,M)

            for ind = 1:M
                ∂ws[ind,ind] = -ns  # diagonal entries

                # Calculate the column index of the off-diagonal entry in the row `ind`.
                sub′ .= ind2sub(N.data, ind)
                if ns == 1  # forward difference
                    if sub′[nw] == Nw
                        sub′[nw] = 1
                    else
                        sub′[nw] += 1
                    end
                else  # ns = -1: backward difference
                    if sub′[nw] == 1
                        sub′[nw] = Nw
                    else
                        sub′[nw] -= 1
                    end
                end

                ind′ = sub2ind(N.data, sub′...)
                ∂ws[ind, ind′] += ns  # off-diagonal entry
            end
            @test create_∂(nw, ns, N) == ∂ws
        end
    end
end  # @testset "create_∂"

N = SVector(3,4,5)
M = prod(N)
r = reshape(collect(1:3M), M, 3)'[:]  # index mapping from block matrix to narrowly banded matrix

@testset "create_curl for U" begin
    # Construct Cu for a uniform grid and BLOCH boundaries.
    ∆ldual = ones.(N.data)
    ebc =  @SVector fill(BLOCH, 3)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    Cu = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine the overall coefficients.
    @test all(any(Cu.≠0, 1))  # no zero columns
    @test all(any(Cu.≠0, 2))  # no zero rows
    @test all(sum(Cu, 2) .== 0)  # all row sums are zero, because Cu * ones(M) = 0

    # Construct Cu for a nonuniform grid and general boundaries.
    ∆ldual = rand.(N.data)
    ebc = SVector(BLOCH, PPC, PDC)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    Cu = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine diagonal blocks.
    @test all(Cu[1:M, 1:M] .== 0)  # (1,1) block
    @test all(Cu[M+1:2M, M+1:2M] .== 0)  # (2,2) block
    @test all(Cu[2M+1:3M, 2M+1:3M] .== 0)  # (3,3) block

    # Examine off-diagonal blocks.
    for nv = nXYZ  # Cartesian component of output vector
        parity = 1
        for nw = next2(nv)  # direction of differentiation
            nw′ = 6 - nv - nw  # Cartesian component of input vector
            ∂w = create_∂(nw, 1, N, ∆ldual[nw], ebc[nw], e⁻ⁱᵏᴸ[nw])

            @test Cu[M*(nv-1)+1:M*nv, M*(nw′-1)+1:M*nw′] == parity * ∂w

            parity = -1
        end
    end

    # Examine reordering.
    Cu_reorder = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=true)
    @test Cu_reorder == Cu[r,r]
end  # @testset "create_curl for U"

@testset "create_curl for V" begin
    # Construct Cv for a uniform grid and BLOCH boundaries.
    ∆lprim = ones.(N.data)
    ebc =  @SVector fill(BLOCH, 3)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    Cv = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine the overall coefficients.
    @test all(any(Cv.≠0, 1))  # no zero columns
    @test all(any(Cv.≠0, 2))  # no zero rows
    @test all(sum(Cv, 2) .== 0)  # all row sums are zero, because Cv * ones(sum(Min)) = 0

    # Construct Cv for a nonuniform grid and general boundaries.
    ∆lprim = rand.(N.data)
    ebc = SVector(BLOCH, PPC, PDC)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    Cv = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine diagonal blocks.
    @test all(Cv[1:M, 1:M] .== 0)  # (1,1) block
    @test all(Cv[M+1:2M, M+1:2M] .== 0)  # (2,2) block
    @test all(Cv[2M+1:3M, 2M+1:3M] .== 0)  # (3,3) block

    # Examine off-diagonal blocks.
    for nv = nXYZ  # Cartesian component of output vector
        parity = 1
        for nw = next2(nv)  # direction of differentiation
            nw′ = 6 - nv - nw  # Cartesian component of input vector
            ∂w = create_∂(nw, -1, N, ∆lprim[nw], ebc[nw], e⁻ⁱᵏᴸ[nw])

            @test Cv[M*(nv-1)+1:M*nv, M*(nw′-1)+1:M*nw′] == parity * ∂w

            parity = -1
        end
    end

    # Examine reordering
    Cv_reorder = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=true)
    @test Cv_reorder == Cv[r,r]
end  # @testset "create_curl for V"

@testset "curl of curl" begin
    # Construct Cu and Cv for a uniform grid and BLOCH boundaries.
    ∆ldual = ones.(N.data)
    ∆lprim = ones.(N.data)
    ebc =  @SVector fill(BLOCH, 3)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    Cu = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)
    Cv = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Test symmetry of each block.
    for i = nXYZ
        for j = next2(i)
            -Cv[(i-1)*M+1:i*M,(j-1)*M+1:j*M]' == Cu[(i-1)*M+1:i*M,(j-1)*M+1:j*M]
        end
    end

    # Construct Cv * Cu
    A = Cv * Cu

    # Test curl of curl.
    @test all(diag(A) .== 4)  # all diagonal entries are 4
    @test all(sum(A.≠0, 2) .== 13)  # 13 nonzero entries per row
    @test A == A'  # Hermitian

    B = A - 4*speye(A)
    @test all(abs.(B[B.≠0]).==1)  # all nonzero off-diagonal entries are ±1
end  # @testset "curl of curl"

end  # @testset "differential"
