@testset "mean" begin

@testset "create_mean" begin
    N = SVector(8,9,10)
    # N = SVector(3,3,3)
    M = prod(N)
    for isfwd = (true, false)
        Mdiag = spzeros(3M,3M)
        Msup = spzeros(3M,3M)
        Msub = spzeros(3M,3M)

        for nw = nXYZ
            Nw = N[nw]
            sub′ = Vector{Int}(undef, 3)

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

            Is = 1+(nw-1)*M:nw*M
            Mdiag[Is,Is] .= Mws

            Js = 1+mod(nw,3)*M:mod1(nw+1,3)*M
            Msup[Is,Js] .= Mws

            Js = 1+mod(nw-2,3)*M:mod1(nw-1,3)*M
            Msub[Is,Js] .= Mws
        end  # for nw
        @test create_mean([isfwd,isfwd,isfwd], N, reorder=false) == Mdiag
        @test create_mean([isfwd,isfwd,isfwd], N, kdiag=1, reorder=false) == Msup
        @test create_mean([isfwd,isfwd,isfwd], N, kdiag=-1, reorder=false) == Msub
    end
end  # @testset "create_m"

# @testset "create_mean" begin
# create_mean(isfwd::AbsVecBool,  # isfwd[w] = true|false for forward|backward averaging
#             N::AbsVecInteger,  # size of grid
#             isbloch::AbsVecBool=fill(true,length(N)),  # boundary conditions in x, y, z
#             e⁻ⁱᵏᴸ::AbsVecNumber=ones(length(N));  # Bloch phase factor in x, y, z
#             kdiag::Integer=0,  # 0|+1|-1 for diagonal|superdiagonal|subdiagonal of material parameter
#             reorder::Bool=true) =  # true for more tightly banded matrix
#
#     N = SVector(8,9,10)
#     # N = SVector(3,3,3)
#     M = prod(N)
#     for nw = nXYZ
#     # for nw = (nX,)
#         Nw = N[nw]
#         sub′ = Vector{Int}(undef, 3)
#
#         for isfwd = (true, false)
#         # for isfwd = (true,)
#             Mws = spzeros(M,M)
#
#             for ind = 1:M
#                 Mws[ind,ind] = 0.5  # diagonal entries
#
#                 # Calculate the column index of the off-diagonal entry in the row `ind`.
#                 sub′ .= CartesianIndices(N.data)[ind].I  # subscripts of off-diagonal entry
#                 if isfwd  # forward difference
#                     if sub′[nw] == Nw
#                         sub′[nw] = 1
#                     else
#                         sub′[nw] += 1
#                     end
#                 else  # isfwd = false (backward difference)
#                     if sub′[nw] == 1
#                         sub′[nw] = Nw
#                     else
#                         sub′[nw] -= 1
#                     end
#                 end
#
#                 ind′ = LinearIndices(N.data)[sub′...]  # index of off-diagonal entry
#                 Mws[ind, ind′] = 0.5  # off-diagonal entry
#             end
#             @test create_m(nw, isfwd, N) == Mws
#         end
#     end
#
# end  # @testset "create_mean"

end  # @testset "mean"
