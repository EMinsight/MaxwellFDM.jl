# The following without the qualifier MaxwellFDM cannot find the Grid constructor of MaxwellFDM.
# Grid(axis::Axis, unit::PhysUnit, lprim::AbsVecReal, Npml::NTuple{2,Int}, isbloch::Bool) =
#     MaxwellFDM.Grid((axis,), unit, (lprim,), ([Npml[nN]], [Npml[nP]]), (isbloch,))

# Calculate ghost points from l, L, and isbloch.
lghost(l::NTuple{2,NTuple{K,AbstractVector{<:Real}}},  # grid point locations
       L::SVector{K,Float64},  # domain size
       isbloch::SVector{K,Bool}  # boundary condition
      ) where {K} =
    (map((lprimₖ,Lₖ) -> lprimₖ[1]+Lₖ, SVector(l[nPR]), L),  # lg[PRIM]
     map((isblochₖ,lprimₖ,ldualₖ,Lₖ) -> (isblochₖ ? ldualₖ[end]-Lₖ : 2lprimₖ[1]-ldualₖ[1]), isbloch, SVector(l[nPR]), SVector(l[nDL]), L))  # lg[DUAL]

@testset "grid" begin
L₀ = 1e-9
unit = PhysUnit(L₀)

@testset "Grid{1}, Bloch boundary" begin
    isbloch = true
    N = 22
    ∆ldual = rand(N)
    L = sum(∆ldual)
    l₀ = L / 2
    lprim = cumsum([-l₀; ∆ldual])

    g1 = Grid(X̂, unit, lprim, isbloch)

    @test g1.unit == unit
    @test g1.N == [N]
    @test g1.L ≈ [L]
    @test all(SVector.(lprim) .∈ Ref(g1))  # `∈` supports only SVector, not Vector
    ldual = MaxwellFDM.movingavg(lprim)
    pop!(lprim)
    @test g1.l ≈ ((lprim,), (ldual,))
    ∆lprim = [ldual[1]+L-ldual[end]; diff(ldual)]
    @test g1.∆l ≈ ((∆lprim,), (∆ldual,))
    @test g1.isbloch == [isbloch]
    # @test g1.Npml == ([Npmln], [Npmlp])
    # @test g1.Lpml ≈ ([sum(∆ldual[1:Npmln])], [sum(∆ldual[end-Npmlp+1:end])])
    # @test g1.lpml ≈ ([sum(∆ldual[1:Npmln])-l₀], [sum(∆ldual)-sum(∆ldual[end-Npmlp+1:end])-l₀])
    @test g1.bounds ≈ ([lprim[1]], [lprim[end]+∆ldual[end]])
    @test SVector(-l₀) ∈ g1  # `∈` supports only SVector
    lg = lghost(((lprim,),(ldual,)), SVector(L), SVector(isbloch))
    lprim_g = g1.ghosted.l[nPR][1]
    ldual_g = g1.ghosted.l[nDL][1]
    @test lprim_g[1:end-1]≈lprim && lprim_g[end]≈lg[nPR][1] && ldual_g[2:end]≈ldual && ldual_g[1]≈lg[nDL][1]
    τlprim_g = g1.ghosted.τl[nPR][1]
    τldual_g = g1.ghosted.τl[nDL][1]
    @test τlprim_g[1:end-1]≈lprim && τlprim_g[end]≈lg[nPR][1]-L && τldual_g[2:end]≈ldual && τldual_g[1]≈lg[nDL][1]+L
    τindprim_g = g1.ghosted.τind[nPR][1]
    τinddual_g = g1.ghosted.τind[nDL][1]
    @test all(τindprim_g[1:end-1].==1:N) && τindprim_g[end]==1 && all(τinddual_g[2:end].==2:N+1) && τinddual_g[1]==N+1
    ∆τprim_g = g1.ghosted.∆τ[nPR][1]
    ∆τdual_g = g1.ghosted.∆τ[nDL][1]
    @test iszero(∆τprim_g[1:end-1]) && ∆τprim_g[end]≈-L && iszero(∆τdual_g[2:end]) && ∆τdual_g[1]≈L
end  # @testset "Grid{1}, primal boundary"


@testset "Grid{3}" begin
    N = [29, 23, 6]
    ∆ldual = ntuple(d->rand(N[d]), numel(Axis))
    isbloch = [false, false, true]

    L = SVector(sum.(∆ldual))  # SVec3Float
    l₀ = L ./ 2  # SVec3Float
    lprim = map((v,s)->v.-s, map(x->[0; cumsum(x)], ∆ldual), (l₀...,))  # tuple of vectors

    g3 = Grid(unit, lprim, isbloch)

    @test g3.unit == unit
    @test g3.N == N
    @test g3.L ≈ L
    ldual = MaxwellFDM.movingavg.(lprim)
    pop!.(lprim)
    @test g3.l ≈ (lprim, ldual)
    ∆lprim = diff.(ldual)
    prepend!(∆lprim[nX], 2(ldual[nX][1]-lprim[nX][1]))
    prepend!(∆lprim[nY], 2(ldual[nY][1]-lprim[nY][1]))
    prepend!(∆lprim[nZ], ldual[nZ][1]+L[nZ]-ldual[nZ][end])
    @test g3.∆l ≈ (∆lprim, ∆ldual)
    @test g3.isbloch == isbloch
    # @test g3.Npml == Npml
    # @test g3.Lpml ≈ (
    #     [sum(∆ldual[d][1:Npml[nN][d]]) for d = nXYZ],
    #     [sum(∆ldual[d][end-Npml[nP][d]+1:end]) for d = nXYZ],
    # )
    # @test g3.lpml ≈ (
    #     [sum(∆ldual[d][1:Npml[nN][d]]) - l₀[d] for d = nXYZ],
    #     [sum(∆ldual[d]) - sum(∆ldual[d][end-Npml[nP][d]+1:end]) - l₀[d] for d = nXYZ]
    # )
    @test g3.bounds ≈ (
        [lprim[d][1] for d = nXYZ],
        [lprim[d][end] + ∆ldual[d][end] for d = nXYZ]
    )
    @test -l₀ ∈ g3
    @test all(g3.bounds .∈ Ref(g3))
    lg = lghost((lprim,ldual), L, g3.isbloch)
    lprim_g = g3.ghosted.l[nPR]
    ldual_g = g3.ghosted.l[nDL]
    @test pop!.(lprim_g)≈lg[nPR].data && popfirst!.(ldual_g)≈lg[nDL].data && lprim_g≈lprim && ldual_g≈ldual
    τlprim_g = g3.ghosted.τl[nPR]
    τldual_g = g3.ghosted.τl[nDL]
    @test pop!.(τlprim_g)≈lg[nPR].data.-(0,0,L[nZ]) && popfirst!.(τldual_g)≈(ldual[nX][1],ldual[nY][1],lg[nDL][nZ]+L[nZ]) && τlprim_g≈lprim && τldual_g≈ldual
    τindprim_g = g3.ghosted.τind[nPR]
    τinddual_g = g3.ghosted.τind[nDL]
    @test pop!.(τindprim_g)==(N[nX]+1,N[nY]+1,1) && popfirst!.(τinddual_g)==(2,2,N[nZ]+1) && τindprim_g==map(m->collect(1:m), tuple(N...)) && τinddual_g==map(m->collect(2:m+1), tuple(N...))
    ∆τprim_g = g3.ghosted.∆τ[nPR]
    ∆τdual_g = g3.ghosted.∆τ[nDL]
    @test pop!.(∆τprim_g)≈(0,0,-L[nZ]) && popfirst!.(∆τdual_g)≈(0,0,L[nZ]) && all(iszero.(∆τprim_g)) && all(iszero.(∆τdual_g))
end  # @testset "Grid{3}"

# @testset "Grid{2}" begin
#     isbloch = (false, false, true)
#     Npml = ((10,9), (5,6), (1,4))
#     N = (10, 12, 1)
#     M = (+).(sum.(Npml), N)
#     ∆ldual = rand.(M)  # tuple of arrays
#     L = sum.(∆ldual)
#     l₀ = (/).(L, 2)
#     lprim = (-).(map(x->[0; cumsum(x)], ∆ldual), l₀)
#
#     g3 = Grid3D(unit, lprim, Npml, isbloch)
#     normal_axis = Ŷ
#     h, v = Int.(next3(normal_axis))
#     g2 = Grid2D(g3, normal_axis)
#
#     lprim = (lprim[h], lprim[v])
#     ∆ldual = (∆ldual[h], ∆ldual[v])
#     M = (M[h], M[v])
#     L = (L[h], L[v])
#     isbloch = (isbloch[h], isbloch[v])
#     Npml = (Npml[h], Npml[v])
#     l₀ = (l₀[h], l₀[v])
#
#     @test g2.normal_axis == Ŷ
#     @test unit_(g2) == unit
#     @test g2.N == M
#     @test g2.L ≈ L
#     @test bounds_(g2) ≈ ((lprim[nHRZ][1],lprim[nHRZ][end]), (lprim[nVRT][1],lprim[nVRT][end]))
#     zdual = (lprim[nHRZ][1:end-1] + lprim[nHRZ][2:end]) / 2
#     xdual = (lprim[nVRT][1:end-1] + lprim[nVRT][2:end]) / 2
#     ldual = (zdual, xdual)
#     lprim = (nw->lprim[nw][1:end-1]).(nPD)
#     @test l_(g2) ≈ (lprim, ldual)
#     ∆lprim = (
#         [ldual[nHRZ][1]+L[nHRZ]-ldual[1][end]; diff(zdual)],
#         [2(ldual[nVRT][1]-lprim[nVRT][1]); diff(xdual)]
#     )
#     @test ∆l_(g2) ≈ (∆lprim, ∆ldual)
#     @test isbloch_(g2) == isbloch
#     @test Npml_(g2) == Npml
#     @test Lpml_(g2) ≈ (
#         (sum(∆ldual[nHRZ][1:Npml[nHRZ][nN]]), sum(∆ldual[nHRZ][end-Npml[nHRZ][nP]+1:end])),
#         (sum(∆ldual[nVRT][1:Npml[nVRT][nN]]), sum(∆ldual[nVRT][end-Npml[nVRT][nP]+1:end]))
#     )
#     @test lpml_(g2) ≈ (
#         (sum(∆ldual[nHRZ][1:Npml[nHRZ][nN]]) - l₀[nHRZ], sum(∆ldual[nHRZ]) - sum(∆ldual[nHRZ][end-Npml[nHRZ][nP]+1:end]) - l₀[nHRZ]),
#         (sum(∆ldual[nVRT][1:Npml[nVRT][nN]]) - l₀[nVRT], sum(∆ldual[nVRT]) - sum(∆ldual[nVRT][end-Npml[nVRT][nP]+1:end]) - l₀[nVRT])
#     )
#     @test contains(g2, map(-,l₀))
#     @test all(.!contains.(g2, [bounds_(g2,NEG), bounds_(g2,POS)], false))
#     @test all(contains.(g2, [bounds_(g2,NEG), bounds_(g2,POS)]))
#     # @test pbounds(g2, withpml=true) ≈ ((-l₀[1],sum(∆ldual[1])-l₀[1]), (-l₀[2],sum(∆ldual[2])-l₀[2]))
#     # @test pbounds(g2, withpml=false) ≈ lpml_(g2)
#     # @test lplot(g2, withinterp=true, withpml=true) ≈ (
#     #     (lprim[1],[lprim[1][1];ldual[1][2:end-1];lprim[1][end]]),
#     #     (lprim[2],[lprim[2][1];ldual[2][2:end-1];lprim[2][end]])
#     # )
#     # @test voxel_edges(g2, withpml=true) ≈ ((ldual[1],lprim[1]), (ldual[2],lprim[2]))
# end  # @testset "Grid{2}"

@testset "isproper_blochphase" begin
    @test isproper_blochphase(1im, true)
    @test isproper_blochphase(cos(π/0.9)+sin(π/0.9)im, true)
    @test !isproper_blochphase(2., true)
    @test !isproper_blochphase(2., false)
end

end  # @testset "grid"
