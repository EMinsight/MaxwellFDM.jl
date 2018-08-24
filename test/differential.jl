@testset "differential" begin

@testset "create_∂" begin
    N = SVector(8,9,10)
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
            @test create_∂(nw, ns==1, N) == ∂ws
        end
    end
end  # @testset "create_∂"

N = SVector(3,4,5)
M = prod(N)
r = reshape(collect(1:3M), M, 3)'[:]  # index mapping from block matrix to narrowly banded matrix
Z = spzeros(M,M)

@testset "create_curl for U" begin
    # Construct Cu for a uniform grid and BLOCH boundaries.
    Cu = create_curl([true,true,true], N, reorder=false)

    # Examine the overall coefficients.
    @test all(any(Cu.≠0, 1))  # no zero columns
    @test all(any(Cu.≠0, 2))  # no zero rows
    @test all(sum(Cu, 2) .== 0)  # all row sums are zero, because Cu * ones(M) = 0

    ∂x = (nw = 1; create_∂(nw, true, N))
    ∂y = (nw = 2; create_∂(nw, true, N))
    ∂z = (nw = 3; create_∂(nw, true, N))
    @test Cu == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Construct Cu for a nonuniform grid and general boundaries.
    ∆ldual = rand.(N.data)
    isbloch = SVector(true, false, false)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    Cu = create_curl([true,true,true], N, ∆ldual, isbloch, e⁻ⁱᵏᴸ, reorder=false)

    # Examine Cu.
    ∂x = (nw = 1; create_∂(nw, true, N, ∆ldual[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw]))
    ∂y = (nw = 2; create_∂(nw, true, N, ∆ldual[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw]))
    ∂z = (nw = 3; create_∂(nw, true, N, ∆ldual[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw]))
    @test Cu == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Examine reordering.
    Cu_reorder = create_curl([true,true,true], N, ∆ldual, isbloch, e⁻ⁱᵏᴸ, reorder=true)
    @test Cu_reorder == Cu[r,r]
end  # @testset "create_curl for U"

@testset "create_curl for V" begin
    # Construct Cv for a uniform grid and BLOCH boundaries.
    Cv = create_curl([false,false,false], N, reorder=false)

    # Examine the overall coefficients.
    @test all(any(Cv.≠0, 1))  # no zero columns
    @test all(any(Cv.≠0, 2))  # no zero rows
    @test all(sum(Cv, 2) .== 0)  # all row sums are zero, because Cv * ones(sum(Min)) = 0

    ∂x = (nw = 1; create_∂(nw, false, N))
    ∂y = (nw = 2; create_∂(nw, false, N))
    ∂z = (nw = 3; create_∂(nw, false, N))
    @test Cv == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Construct Cv for a nonuniform grid and general boundaries.
    ∆lprim = rand.(N.data)
    isbloch = SVector(true, false, false)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    Cv = create_curl([false,false,false], N, ∆lprim, isbloch, e⁻ⁱᵏᴸ, reorder=false)

    # Examine Cv.
    ∂x = (nw = 1; create_∂(nw, false, N, ∆lprim[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw]))
    ∂y = (nw = 2; create_∂(nw, false, N, ∆lprim[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw]))
    ∂z = (nw = 3; create_∂(nw, false, N, ∆lprim[nw], isbloch[nw], e⁻ⁱᵏᴸ[nw]))
    @test Cv == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Examine reordering
    Cv_reorder = create_curl([false,false,false], N, ∆lprim, isbloch, e⁻ⁱᵏᴸ, reorder=true)
    @test Cv_reorder == Cv[r,r]
end  # @testset "create_curl for V"

@testset "curl of curl" begin
    # Construct Cu and Cv for a uniform grid and BLOCH boundaries.
    ∆ldual = ones.(N.data)
    ∆lprim = ones.(N.data)
    isbloch =  SVector(true, false, false)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    Cu = create_curl([true,true,true], N, ∆ldual, isbloch, e⁻ⁱᵏᴸ, reorder=false)
    Cv = create_curl([false,false,false], N, ∆lprim, isbloch, e⁻ⁱᵏᴸ, reorder=false)

    # Test symmetry of each block.
    for i = nXYZ
        for j = next2(i)
            -Cv[(i-1)*M+1:i*M,(j-1)*M+1:j*M]' == Cu[(i-1)*M+1:i*M,(j-1)*M+1:j*M]
        end
    end

    # Construct Cv * Cu for all BLOCH.
    isbloch =  @SVector fill(true, 3)
    Cu = create_curl([true,true,true], N, ∆ldual, isbloch, e⁻ⁱᵏᴸ, reorder=false)
    Cv = create_curl([false,false,false], N, ∆lprim, isbloch, e⁻ⁱᵏᴸ, reorder=false)
    A = Cv * Cu

    # Test curl of curl.
    @test all(diag(A) .== 4)  # all diagonal entries are 4
    @test all(sum(A.≠0, 2) .== 13)  # 13 nonzero entries per row
    @test A == A'  # Hermitian

    B = A - 4*speye(A)
    @test all(abs.(B[B.≠0]).==1)  # all nonzero off-diagonal entries are ±1
end  # @testset "curl of curl"

@testset "curl of curl, mixed forward and backward" begin
    # Construct Cu and Cv for a uniform grid and BLOCH boundaries.
    ∆ldual = ones.(N.data)
    ∆lprim = ones.(N.data)
    isbloch =  SVector(true, false, false)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    isfwd = [true,false,true]
    Cu = create_curl(isfwd, N, ∆ldual, isbloch, e⁻ⁱᵏᴸ, reorder=false)
    Cv = create_curl(.!isfwd, N, ∆lprim, isbloch, e⁻ⁱᵏᴸ, reorder=false)

    # Test symmetry of each block.
    for i = nXYZ
        for j = next2(i)
            -Cv[(i-1)*M+1:i*M,(j-1)*M+1:j*M]' == Cu[(i-1)*M+1:i*M,(j-1)*M+1:j*M]
        end
    end

    # Construct Cv * Cu for all BLOCH.
    isbloch =  @SVector fill(true, 3)
    Cu = create_curl(isfwd, N, ∆ldual, isbloch, e⁻ⁱᵏᴸ, reorder=false)
    Cv = create_curl(.!isfwd, N, ∆lprim, isbloch, e⁻ⁱᵏᴸ, reorder=false)
    A = Cv * Cu

    # Test curl of curl.
    @test all(diag(A) .== 4)  # all diagonal entries are 4
    @test all(sum(A.≠0, 2) .== 13)  # 13 nonzero entries per row
    @test A == A'  # Hermitian

    B = A - 4*speye(A)
    @test all(abs.(B[B.≠0]).==1)  # all nonzero off-diagonal entries are ±1
end  # @testset "curl of curl"

end  # @testset "differential"
