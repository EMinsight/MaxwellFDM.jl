@testset "pml" begin

@testset "PMLParam" begin
    pml = PMLParam()
    @test pml.m===4.0 && pml.R===exp(-16) && pml.κmax===1.0 && pml.amax===0.0 && pml.ma===4.0
end  # @testset "PMLParam"

@testset "1D" begin
    Npmln = 5
    Npmlp = 7
    Npml = (Npmln, Npmlp)
    N = 10
    M = N + Npmln + Npmlp
    ∆ldual = rand(M)
    L = sum(∆ldual)
    l₀ = L / 2
    lprim = cumsum([-l₀; ∆ldual])
    ldual = MaxwellFDM.movingavg(lprim)

    lpml = (lprim[1+Npmln], lprim[end-Npmlp])
    Lpml = (lpml[nN]-lprim[1], lprim[end]-lpml[nP])

    @test get_pml_loc(lprim, Npml) == (lpml, Lpml)

    ω = rand()
    pml = PMLParam()
    σmax = -(pml.m+1) * log(pml.R) / 2 ./ Lpml  # see calc_stretch_factor
    f(l) = l≤lpml[nN] ? (1 + σmax[nN] * ((lpml[nN]-l)/Lpml[nN])^pml.m / (im*ω)) : (l≥lpml[nP] ? (1 + σmax[nP] * ((l-lpml[nP])/Lpml[nP])^pml.m / (im*ω)) : 1)
    sprim = f.(lprim[1:end-1])
    sdual = f.(ldual)

    @test gen_stretch_factor(ω, (lprim[1:end-1], ldual), lpml, Lpml) == (sprim, sdual)
end  # @testset "get_pml_loc"

@testset "1D, Npml = 0" begin
    Npml = (0,0)
    lprim = [0.0, 1.0]  # with ghost point
    ldual = [0.5]  # without ghost point

    lpml = (0.0, 1.0)
    Lpml = (0.0, 0.0)

    @test get_pml_loc(lprim, Npml) == (lpml, Lpml)

    ω = rand()
    @test gen_stretch_factor(ω, (lprim[1:end-1], ldual), lpml, Lpml) == ([1.0], [1.0])
end  # @testset "get_pml_loc"

@testset "3D" begin
    Npml = ([10,5,1], [9,6,4])
    N = [10, 12, 1]
    M = sum(Npml) + N
    ∆ldual = ntuple(d->rand(M[d]), numel(Axis))

    L = SVector(sum.(∆ldual))  # SVec3Float
    l₀ = L ./ 2  # SVec3Float
    lprim = map((v,s)->v.-s, map(x->[0; cumsum(x)], ∆ldual), (l₀...,))  # tuple of vectors
    ldual = MaxwellFDM.movingavg.(lprim)

    lpml = (SVector{3}((w->lprim[w][1+Npml[nN][w]]).(nXYZ)), SVector{3}((w->lprim[w][end-Npml[nP][w]]).(nXYZ)))
    Lpml = (lpml[nN].-SVector{3}((w->lprim[w][1]).(nXYZ)), SVector{3}((w->lprim[w][end]).(nXYZ)).-lpml[nP])

    @test get_pml_loc(lprim, Npml) == (lpml, Lpml)

    ω = rand()
    pml = PMLParam()
    σmax = (-(pml.m+1) * log(pml.R) / 2 ./ Lpml[nN], -(pml.m+1) * log(pml.R) / 2 ./ Lpml[nP])  # see calc_stretch_factor
    f(l, nw) = l≤lpml[nN][nw] ? (1 + σmax[nN][nw] * ((lpml[nN][nw]-l)/Lpml[nN][nw])^pml.m / (im*ω)) : (l≥lpml[nP][nw] ? (1 + σmax[nP][nw] * ((l-lpml[nP][nw])/Lpml[nP][nw])^pml.m / (im*ω)) : 1)
    sprim = (f.(lprim[nX][1:end-1],nX), f.(lprim[nY][1:end-1],nY), f.(lprim[nZ][1:end-1],nZ))
    sdual = (f.(ldual[nX],nX), f.(ldual[nY],nY), f.(ldual[nZ],nZ))

    @test gen_stretch_factor(ω, ((l->l[1:end-1]).(lprim), ldual), lpml, Lpml) == (sprim, sdual)
end

@testset "3D, Npml = 0" begin
    Npml = ([0,0,0], [0,0,0])
    lprim = ([0.0,1.0], [0.0,1.0], [0.0,1.0])  # with ghost points
    ldual = ([0.5], [0.5], [0.5])  # without ghost points

    lpml = ([0.0,0.0,0.0], [1.0,1.0,1.0])
    Lpml = ([0.0,0.0,0.0], [0.0,0.0,0.0])

    @test get_pml_loc(lprim, Npml) == (lpml, Lpml)

    ω = rand()
    @test gen_stretch_factor(ω, ((l->l[1:end-1]).(lprim), ldual), lpml, Lpml) == (([1.0],[1.0],[1.0]), ([1.0],[1.0],[1.0]))
end  # @testset "get_pml_loc"

end  # @testset "pml"
