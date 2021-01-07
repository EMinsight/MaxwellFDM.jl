@testset "source" begin

nX, nY, nZ = 1, 2, 3

@testset "distweights" begin
    ∆lg = 2:2:22  # make sure ∆l's are varying to prevent false positive
    l′ = cumsum([-1; ∆lg])  # [-1,1,5,11,19,29,41,55,71,89,109,131]: spacing is even

    lprim_g = StaggeredGridCalculus.movingavg(l′)  # [0,3,8,15,24,35,48,63,80,99,120]: spacing is odd, 11 entries (l has 10 entry)
    ldual_g = l′[1:end-1]  # [-1,1,5,11,19,29,41,55,71,89,109]: spacing is even

    domain = [lprim_g[1], lprim_g[end]]

    lprim, ∆lprim = lprim_g[1:end-1], diff(ldual_g)
    ldual, ∆ldual = ldual_g[2:end], diff(lprim_g)

    # Test for sources to assign at the primal grid points [0,3,8,15,24,35,48,63,80,99 (,120)].
    c = 8  # at grid point
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([3,3], [1/∆lprim[3], 0]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([3,3], [1/∆lprim[3], 0]))

    c = 11  # between grid points
    r = (c-8) / (15-8)
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([3,4], [1/∆lprim[3]*(1-r), 1/∆lprim[4]*r]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([3,4], [1/∆lprim[3]*(1-r), 1/∆lprim[4]*r]))

    c = 3  # farthest point affected by negative-end boundary
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([2,2], [1/∆lprim[2], 0]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([2,2], [1/∆lprim[2], 0]))

    c = 1  # within range affected by negative-end boundary
    r = (3-c) / (3-0)
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([1,2], [1/∆lprim[1]*r, 1/∆lprim[2]*(1-r)]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([2,2], [1/∆lprim[2]*(1-r), 0]))

    c = 0  # negative-end boundary
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([1,1], [1/∆lprim[1], 0]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([2,2], [0, 0]))

    c = -0.5  # within range affected by negative-end PDC, but outside Bloch and PPC domain (starting at l = 0)
    @test_throws ArgumentError MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true)
    @test_throws ArgumentError MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false)

    c = -1  # behind negative-end boundary
    @test_throws ArgumentError MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true)
    @test_throws ArgumentError MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false)

    c = 99  # farthest point affected by positive-end boundary
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([10,10], [1/∆lprim[10], 0]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([10,10], [1/∆lprim[10], 0]))

    c = 105  # within range affected by positive-end boundary
    r = (c-99) / (120-99)
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([10,1], [1/∆lprim[10]*(1-r), 1/∆lprim[1]*r]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([10,10], [1/∆lprim[10]*(1-r), 0]))

    c = 120  # positive-end boundary
    r = (c-99) / (120-99)
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, true) .≈ ([10,1], [1/∆lprim[10]*(1-r), 1/∆lprim[1]*r]))
    @test all(MaxwellFDM.distweights(c, PRIM, domain, lprim, ∆lprim, false) .≈ ([10,10], [1/∆lprim[10]*(1-r), 0]))


    # Test for sources to assign at the dual grid points [(-1,) 1,5,11,19,29,41,55,71,89,109].
    # Note the domain boundaries are [0, 120].
    c = 11  # at grid point
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([3,3], [1/∆ldual[3], 0]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([3,3], [1/∆ldual[3], 0]))

    c = 15  # between grid points
    r = (c-11) / (19-11)
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([3,4], [1/∆ldual[3]*(1-r), 1/∆ldual[4]*r]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([3,4], [1/∆ldual[3]*(1-r), 1/∆ldual[4]*r]))

    c = 1  # farthest point affected by negative-end boundary
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([1,1], [1/∆ldual[1], 0]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([1,1], [1/∆ldual[1], 0]))

    c = 0.5  # within range affected by negative-end boundary
    r = (1-c) / ((1-0) + (120-109))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([1,10], [1/∆ldual[1]*(1-r), 1/∆ldual[10]*r]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([1,1], [1/∆ldual[1], 0]))

    c = 0  # negative-end boundary
    r = (1-c) / ((1-0) + (120-109))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([1,10], [1/∆ldual[1]*(1-r), 1/∆ldual[10]*r]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([1,1], [1/∆ldual[1], 0]))

    c = -0.5  # behind negative-end boundary
    @test_throws ArgumentError MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true)
    @test_throws ArgumentError MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false)

    c = 109  # farthest point affected by positive-end boundary
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([10,10], [1/∆ldual[10], 0]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([10,10], [1/∆ldual[10], 0]))

    c = 115  # within range affected by positive-end boundary
    r = (c-109) / ((120-109) + (1-0))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([10,1], [1/∆ldual[10]*(1-r), 1/∆ldual[1]*r]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([10,10], [1/∆ldual[10], 0]))

    c = 120  # positive-end boundary
    r = (c-109) / ((120-109) + (1-0))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, true) .≈ ([10,1], [1/∆ldual[10]*(1-r), 1/∆ldual[1]*r]))
    @test all(MaxwellFDM.distweights(c, DUAL, domain, ldual, ∆ldual, false) .≈ ([10,10], [1/∆ldual[10], 0]))

    # Test for N = 1
    domain = [0, 2]
    lprim, ldual = [0], [1]
    ∆lprim, ∆ldual = [9], [10]  # this doesn't make sense geometrically, but for test purpose (they are used to get current density from current)
    c₁, c₂ = 2rand(), 2rand()  # within domain

    ind_c₁, wt_c₁ = MaxwellFDM.distweights(c₁, PRIM, domain, lprim, ∆lprim, true)
    ind_c₂, wt_c₂ = MaxwellFDM.distweights(c₂, PRIM, domain, lprim, ∆lprim, true)
    @test ind_c₁ == ind_c₂ == [1,1]
    @test sum(wt_c₁) ≈ sum(wt_c₂) ≈ 1 / ∆lprim[1]

    ind_c₁, wt_c₁ = MaxwellFDM.distweights(c₁, DUAL, domain, ldual, ∆ldual, true)
    ind_c₂, wt_c₂ = MaxwellFDM.distweights(c₂, DUAL, domain, ldual, ∆ldual, true)
    @test ind_c₁ == ind_c₂ == [1,1]
    @test sum(wt_c₁) ≈ sum(wt_c₂) ≈ 1 / ∆ldual[1]
end  # @testset "distweights"

@testset "PlaneSrc" begin
    src = PlaneSrc([0,0,1], 0, [1,0,0])  # x-polarized source in z-normal plane

    isbloch = [true, true, true]

    # Coarse grid
    lprim = (-10:10, -10:10, -10:10)
    g3 = Grid(lprim, isbloch)
    ∆a = 1.0^2  # area element

    ft = EE
    boundft = SVector(EE,EE,EE)

    j3d = create_field_array(g3.N)
    add!(j3d, ft2gt.(ft,boundft), g3.bounds, g3.l, g3.∆l, g3.isbloch, src)

    @test maximum(abs, j3d) == 1.0  # J = K / ∆z
    @test all(j3d[:,:,:,nY] .== 0)  # Jy = 0
    @test all(j3d[:,:,:,nZ] .== 0)  # Jz = 0

    # Fine grid
    lprim = (-10:0.5:10, -10:0.5:10, -10:0.5:10)
    g3_fine = Grid(lprim, isbloch)

    j3d_fine = create_field_array(g3_fine.N)
    add!(j3d_fine, ft2gt.(ft,boundft), g3_fine.bounds, g3_fine.l, g3_fine.∆l, g3_fine.isbloch, src)
    ∆a_fine = 0.5^2  # area element

    @test sum(j3d[1,:,:,nX]) * ∆a ≈ sum(j3d_fine[1,:,:,nX]) * ∆a_fine  # total current through x-normal cross section is independent of grid resolution
    @test maximum(abs, j3d_fine) == 2.0  # J = K / ∆z

    # Nonuniform grid in z
    zprim = sort(rand(21)) * 20
    z_avg = mean(zprim)
    zprim .-= z_avg  # z = 0 is within the z-range
    lprim = (-10:10, -10:10, zprim)
    g3_nu = Grid(lprim, isbloch)

    ∆yprim = g3_nu.∆l[nPR][nY]
    ∆zprim = g3_nu.∆l[nPR][nZ]

    j3d_nu = create_field_array(g3_nu.N)
    add!(j3d_nu, ft2gt.(ft,boundft), g3_nu.bounds, g3_nu.l, g3_nu.∆l, g3_nu.isbloch, src)

    @test sum(j3d[1,:,:,nX]) * ∆a ≈ sum(j3d_nu[1,:,:,nX] .* (∆yprim * transpose(∆zprim)))  # total current through x-normal cross section is independent of grid resolution

    # Verify if the total current does not change with the location of the intercept in 3D for Bloch.
end  # @testset "PlaneSrc"

@testset "PointSrc" begin
    src = PointSrc([0.7,0.7,0.7], [1,1,1])  # polarized in [1,1,1] direction

    isbloch = [true, true, true]

    # Coarse grid
    lprim = (-10:10, -10:10, -10:10)
    g3 = Grid(lprim, isbloch)
    ∆v = 1.0^3  # volume element (dipole-normal area element * dipole length)

    ft = EE
    boundft = SVector(EE,EE,EE)

    j3d = create_field_array(g3.N)
    add!(j3d, ft2gt.(ft,boundft), g3.bounds, g3.l, g3.∆l, g3.isbloch, src)

    @test sum(j3d[:,:,:,nX] .!= 0) == sum(j3d[:,:,:,nY] .!= 0) == sum(j3d[:,:,:,nZ] .!= 0) == 8
    @test sum(j3d[:,:,:,nX]) * ∆v == src.I∆r * src.p[nX]
    @test sum(j3d[:,:,:,nY]) * ∆v == src.I∆r * src.p[nY]
    @test sum(j3d[:,:,:,nZ]) * ∆v == src.I∆r * src.p[nZ]

    # Fine grid
    lprim = (-10:0.5:10, -10:0.5:10, -10:0.5:10)
    g3_fine = Grid(lprim, isbloch)

    j3d_fine = create_field_array(g3_fine.N)
    add!(j3d_fine, ft2gt.(ft,boundft), g3_fine.bounds, g3_fine.l, g3_fine.∆l, g3_fine.isbloch, src)
    ∆v_fine = 0.5^3  # volume element (dipole-normal area element * dipole length)

    @test sum(j3d_fine[:,:,:,nX] .!= 0) == sum(j3d_fine[:,:,:,nY] .!= 0) == sum(j3d_fine[:,:,:,nZ] .!= 0) == 8
    @test sum(j3d_fine[:,:,:,nX]) * ∆v_fine == src.I∆r * src.p[nX]
    @test sum(j3d_fine[:,:,:,nY]) * ∆v_fine == src.I∆r * src.p[nY]
    @test sum(j3d_fine[:,:,:,nZ]) * ∆v_fine == src.I∆r * src.p[nZ]


    # Nonuniform grid in x, y, z
    xprim = sort(rand(21)) * 20; x_avg = mean(xprim); xprim .-= (x_avg-0.7)  # x = 0.7 is within the x-range (average x location is moved to 0.7)
    yprim = sort(rand(21)) * 20; y_avg = mean(yprim); yprim .-= (y_avg-0.7)  # y = 0.7 is within the y-range (average y location is moved to 0.7)
    zprim = sort(rand(21)) * 20; z_avg = mean(zprim); zprim .-= (z_avg-0.7)  # z = 0.7 is within the z-range (average z location is moved to 0.7)
    lprim = (xprim, yprim, zprim)
    g3_nu = Grid(lprim, isbloch)

    ∆xprim = g3_nu.∆l[nPR][nX]
    ∆yprim = g3_nu.∆l[nPR][nY]
    ∆zdual = g3_nu.∆l[nDL][nZ]

    j3d_nu = create_field_array(g3_nu.N)
    add!(j3d_nu, ft2gt.(ft,boundft), g3_nu.bounds, g3_nu.l, g3_nu.∆l, g3_nu.isbloch, src)

    Nx, Ny, Nz = g3_nu.N.data
    @test sum(j3d[:,:,:,nZ]) * ∆v ≈ sum(j3d_nu[:,:,:,nZ] .* (reshape(∆xprim,Nx,1,1) .* reshape(∆yprim,1,Ny,1) .* reshape(∆zdual,1,1,Nz)))  # total current through x-normal cross section is independent of grid resolution

    # Verify if the total current dipole strength does not change with the location of the point
    # source in 3D for Bloch.

    # Verify if the total current dipole strength does not change with the z-location in 2D
    # assignment for Bloch (because the current assigned at the positive boundary is added to the
    # negative boundary)
end  # @testset "PlaneSrc"

end  # @testset "source"
