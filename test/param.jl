@testset "param" begin

# Need to test symmetry of param3dmat for different boundary conditions.
@testset "symmetry of create_paramop, uniform grid" begin
    # Create a grid.
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5])
    g3 = Grid(lprim, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = Material{3,3}("Vacuum", ε=εvac)

    εdiel = 2.0
    diel = Material{3,3}("Dielectric", ε=εdiel)

    # Create objects.
    dom_vac = Object(Box(g3.bounds), vac)
    obj_diel = Object(Box([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Sphere([0,0,0], 1), diel)

    # Add objects.
    oind2shp = Shape3[]
    oind2εind = ParamInd[]
    oind2μind = ParamInd[]
    εind2ε = SSComplex3[]
    μind2μ = SSComplex3[]

    add!(oind2shp, (oind2εind,oind2μind), (εind2ε,μind2μ), dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    N = g3.N
    ε3d = create_param_array(N)
    εxx_oind3d = create_oind_array(N)
    εyy_oind3d = create_oind_array(N)
    εzz_oind3d = create_oind_array(N)
    εoo_oind3d = create_oind_array(N)

    μ3d = create_param_array(N)
    μxx_oind3d = create_oind_array(N)
    μyy_oind3d = create_oind_array(N)
    μzz_oind3d = create_oind_array(N)
    μoo_oind3d = create_oind_array(N)

    boundft = SVector(EE,EE,EE)
    assign_param!(ε3d, (μxx_oind3d,μyy_oind3d,μzz_oind3d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g3.ghosted.τl, g3.isbloch)
    assign_param!(ε3d, tuple(μoo_oind3d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g3.ghosted.τl, g3.isbloch)
    assign_param!(μ3d, (εxx_oind3d,εyy_oind3d,εzz_oind3d), ft2gt.(HH,boundft), oind2shp, oind2μind, μind2μ, g3.ghosted.τl, g3.isbloch)
    assign_param!(μ3d, tuple(εoo_oind3d), ft2gt.(HH,boundft), oind2shp, oind2μind, μind2μ, g3.ghosted.τl, g3.isbloch)

    # Perform smoothing.
    smooth_param!(ε3d, (εxx_oind3d,εyy_oind3d,εzz_oind3d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)
    smooth_param!(ε3d, tuple(εoo_oind3d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    ∆l = g3.∆l[nDL]
    ∆l′ = g3.∆l[nPR]
    ∆l′⁻¹ = map(x->inv.(x), ∆l′)

    Mε = create_paramop(ε3d, ft2gt.(EE,boundft), N, ∆l, ∆l′⁻¹, g3.isbloch, order_cmpfirst=true)

    # @test issymmetric(Mε)
    @test Mε ≈ transpose(Mε)

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mε3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "create_paramop"

@testset "symmetry of create_paramop, nonuniform grid" begin
    # Create a grid.
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5])
    g3 = Grid(lprim, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = Material{3,3}("Vacuum", ε=εvac)

    εdiel = 2.0
    diel = Material{3,3}("Dielectric", ε=εdiel)

    # Create objects.
    dom_vac = Object(Box(g3.bounds), vac)
    obj_diel = Object(Box([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Sphere([0,0,0], 1), diel)

    # Add objects.
    oind2shp = Shape3[]
    oind2εind = ParamInd[]
    oind2μind = ParamInd[]
    εind2ε = SSComplex3[]
    μind2μ = SSComplex3[]

    add!(oind2shp, (oind2εind,oind2μind), (εind2ε,μind2μ), dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    N = g3.N
    ε3d = create_param_array(N)
    εxx_oind3d = create_oind_array(N)
    εyy_oind3d = create_oind_array(N)
    εzz_oind3d = create_oind_array(N)
    εoo_oind3d = create_oind_array(N)

    μ3d = create_param_array(N)
    μxx_oind3d = create_oind_array(N)
    μyy_oind3d = create_oind_array(N)
    μzz_oind3d = create_oind_array(N)
    μoo_oind3d = create_oind_array(N)

    boundft = SVector(EE,EE,EE)
    assign_param!(ε3d, (μxx_oind3d,μyy_oind3d,μzz_oind3d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g3.ghosted.τl, g3.isbloch)
    assign_param!(ε3d, tuple(μoo_oind3d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g3.ghosted.τl, g3.isbloch)
    assign_param!(μ3d, (εxx_oind3d,εyy_oind3d,εzz_oind3d), ft2gt.(HH,boundft), oind2shp, oind2μind, μind2μ, g3.ghosted.τl, g3.isbloch)
    assign_param!(μ3d, tuple(εoo_oind3d), ft2gt.(HH,boundft), oind2shp, oind2μind, μind2μ, g3.ghosted.τl, g3.isbloch)

    # Perform smoothing.
    smooth_param!(ε3d, (εxx_oind3d,εyy_oind3d,εzz_oind3d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)
    smooth_param!(ε3d, tuple(εoo_oind3d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    ∆l = g3.∆l[nDL]
    ∆l′ = g3.∆l[nPR]
    ∆l′⁻¹ = map(x->inv.(x), ∆l′)

    Mε = create_paramop(ε3d, ft2gt.(EE,boundft), N, ∆l, ∆l′⁻¹, g3.isbloch, order_cmpfirst=true)

    # Ml
    ∆lprim = g3.∆l[nPR]
    ∆ldual = g3.∆l[nDL]

    Nx, Ny, Nz = N
    nX, nY, nZ = 1, 2, 3
    ∆X = repeat(∆ldual[nX], outer=(1,Ny,Nz))
    ∆Y = repeat(reshape(∆ldual[nY], (1,Ny,1)), outer=(Nx,1,Nz))
    ∆Z = repeat(reshape(∆ldual[nZ], (1,1,Nz)), outer=(Nx,Ny,1))

    ∆L = [reshape(∆X,1,:); reshape(∆Y,1,:); reshape(∆Z,1,:)]
    ∆l = ∆L[:]
    Ml⁻¹ = sparse(Diagonal(1 ./ ∆l))

    # Ma
    ∆X = repeat(∆lprim[nX], outer=(1,Ny,Nz))
    ∆Y = repeat(reshape(∆lprim[nY], (1,Ny,1)), outer=(Nx,1,Nz))
    ∆Z = repeat(reshape(∆lprim[nZ], (1,1,Nz)), outer=(Nx,Ny,1))

    ∆YZ = ∆Y .* ∆Z
    ∆ZX = ∆Z .* ∆X
    ∆XY = ∆X .* ∆Y

    ∆A = [reshape(∆YZ,1,:); reshape(∆ZX,1,:); reshape(∆XY,1,:)]
    ∆a = ∆A[:]
    Ma = sparse(Diagonal(∆a))

    Mε_sym = Ma * Mε * Ml⁻¹

    @test Mε_sym ≈ transpose(Mε_sym)

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mε3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "create_paramop"

# Need to test if PML still works after the field-averaging operator applied to the material
# parameters, because effectively it makes the material inhomogeneous along the direction
# normal to the PML interface.

end  # @testset "param"
