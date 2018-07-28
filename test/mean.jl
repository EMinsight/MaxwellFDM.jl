@testset "mean" begin

@testset "create_m" begin
    N = SVector(8,9,10)
    # N = SVector(3,3,3)
    M = prod(N)
    for nw = nXYZ
    # for nw = (1,)
        Nw = N[nw]
        sub′ = Vector{Int}(3)

        for ns = (-1,1)
        # for ns = (1,)
            Mws = spzeros(M,M)

            for ind = 1:M
                Mws[ind,ind] = 0.5  # diagonal entries

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
                Mws[ind, ind′] = 0.5  # off-diagonal entry
            end
            @test create_m(PRIM, nw, ns, N) == create_m(DUAL, nw, ns, N) == Mws
        end
    end
end  # @testset "create_m"

end  # @testset "mean"
