@testset "mean" begin

@testset "create_m" begin
    N = SVector(8,9,10)
    # N = SVector(3,3,3)
    M = prod(N)
    for nw = nXYZ
    # for nw = (nX,)
        Nw = N[nw]
        sub′ = Vector{Int}(undef, 3)

        for isfwd = (true, false)
        # for isfwd = (true,)
            Mws = spzeros(M,M)

            for ind = 1:M
                Mws[ind,ind] = 0.5  # diagonal entries

                # Calculate the column index of the off-diagonal entry in the row `ind`.
                sub′ .= CartesianIndices(N.data)[ind].I  # subscripts of off-diagonal entry
                if isfwd  # forward difference
                    if sub′[nw] == Nw
                        sub′[nw] = 1
                    else
                        sub′[nw] += 1
                    end
                else  # isfwd = false (backward difference)
                    if sub′[nw] == 1
                        sub′[nw] = Nw
                    else
                        sub′[nw] -= 1
                    end
                end

                ind′ = LinearIndices(N.data)[sub′...]  # index of off-diagonal entry
                Mws[ind, ind′] = 0.5  # off-diagonal entry
            end
            @test create_m(nw, isfwd, N) == Mws
        end
    end
end  # @testset "create_m"

end  # @testset "mean"
