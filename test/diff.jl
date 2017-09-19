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

@testset "create_curlU" begin
    # Construct curlU for a uniform grid and BLOCH boundaries.
    ∆ldual = ones.(N.data)
    ebc =  @SVector fill(BLOCH, 3)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    curlU = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine the overall coefficients.
    @test all(any(curlU.≠0, 1))  # no zero columns
    @test all(any(curlU.≠0, 2))  # no zero rows
    @test all(sum(curlU, 2) .== 0)  # all row sums are zero, because curlU * ones(M) = 0

    # Construct curlU for a nonuniform grid and general boundaries.
    ∆ldual = rand.(N.data)
    ebc = SVector(BLOCH, PPC, PDC)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    curlU = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine diagonal blocks.
    @test all(curlU[1:M, 1:M] .== 0)  # (1,1) block
    @test all(curlU[M+1:2M, M+1:2M] .== 0)  # (2,2) block
    @test all(curlU[2M+1:3M, 2M+1:3M] .== 0)  # (3,3) block

    # Examine off-diagonal blocks.
    for nv = nXYZ  # Cartesian component of output vector
        parity = 1
        for nw = next2(nv)  # direction of differentiation
            nw′ = 6 - nv - nw  # Cartesian component of input vector
            ∂w = create_∂(nw, 1, N, ∆ldual[nw], ebc[nw], e⁻ⁱᵏᴸ[nw])

            @test curlU[M*(nv-1)+1:M*nv, M*(nw′-1)+1:M*nw′] == parity * ∂w

            parity = -1
        end
    end

    # Examine reordering.
    curlU_reorder = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=true)
    @test curlU_reorder == curlU[r,r]
end  # @testset "create_curlU"

@testset "create_curlV" begin
    # Construct curlV for a uniform grid and BLOCH boundaries.
    ∆lprim = ones.(N.data)
    ebc =  @SVector fill(BLOCH, 3)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    curlV = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine the overall coefficients.
    @test all(any(curlV.≠0, 1))  # no zero columns
    @test all(any(curlV.≠0, 2))  # no zero rows
    @test all(sum(curlV, 2) .== 0)  # all row sums are zero, because curlU * ones(sum(Min)) = 0

    # Construct curlU for a nonuniform grid and general boundaries.
    ∆lprim = rand.(N.data)
    ebc = SVector(BLOCH, PPC, PDC)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    curlV = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine diagonal blocks.
    @test all(curlV[1:M, 1:M] .== 0)  # (1,1) block
    @test all(curlV[M+1:2M, M+1:2M] .== 0)  # (2,2) block
    @test all(curlV[2M+1:3M, 2M+1:3M] .== 0)  # (3,3) block

    # Examine off-diagonal blocks.
    for nv = nXYZ  # Cartesian component of output vector
        parity = 1
        for nw = next2(nv)  # direction of differentiation
            nw′ = 6 - nv - nw  # Cartesian component of input vector
            ∂w = create_∂(nw, -1, N, ∆lprim[nw], ebc[nw], e⁻ⁱᵏᴸ[nw])

            @test curlV[M*(nv-1)+1:M*nv, M*(nw′-1)+1:M*nw′] == parity * ∂w

            parity = -1
        end
    end

    # Examine reordering
    curlV_reorder = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=true)
    @test curlV_reorder == curlV[r,r]
end  # @testset "create_curlV"

@testset "curl of curl" begin
    # Construct curlU and curlV for a uniform grid and BLOCH boundaries.
    ∆ldual = ones.(N.data)
    ∆lprim = ones.(N.data)
    ebc =  @SVector fill(BLOCH, 3)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    curlU = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)
    curlV = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Test symmetry of each block.
    for i = nXYZ
        for j = next2(i)
            -curlV[(i-1)*M+1:i*M,(j-1)*M+1:j*M]' == curlU[(i-1)*M+1:i*M,(j-1)*M+1:j*M]
        end
    end

    # Construct curlV * curlU
    A = curlV * curlU

    # Test curl of curl.
    @test all(diag(A) .== 4)  # all diagonal entries are 4
    @test all(sum(A.≠0, 2) .== 13)  # 13 nonzero entries per row
    @test A == A'  # Hermitian

    B = A - 4*speye(A)
    @test all(abs.(B[B.≠0]).==1)  # all nonzero off-diagonal entries are ±1
end  # @testset "curl of curl"

end  # @testset "differential"
