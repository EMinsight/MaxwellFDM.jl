@testset "param" begin

# Need to test symmetry of param3dmat for different boundary conditions.
@testset "symmetry of param3d2mat, uniform grid" begin
    # Create a grid.
    L₀ = 1e-9
    unit = PhysUnit(L₀)
    Npml = ([0,0,0], [0,0,0])
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5])
    g3 = Grid(unit, lprim, Npml, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = EncodedMaterial(PRIM, Material("Vacuum", ε=εvac))

    εdiel = 2.0
    diel = EncodedMaterial(PRIM, Material("Dielectric", ε=εdiel))

    # Create objects.
    dom_vac = Object(Box(g3.bounds), vac)
    obj_diel = Object(Box([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Sphere([0,0,0], 1), diel)

    # Add objects.
    ovec = Object3[]
    paramset = (SMat3Complex[], SMat3Complex[])
    add!(ovec, paramset, dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    param3d = create_param3d(N)
    obj3d = create_n3d(Object3, N)
    pind3d = create_n3d(ParamInd, N)
    oind3d = create_n3d(ObjInd, N)

    assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.isbloch)
    # Test the sanity the assigned param3d here.  It is relatively easy, and it was very helpful.

    smooth_param!(param3d, obj3d, pind3d, oind3d, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    Mξ = param3d2mat(param3d[nPR], PRIM, N, g3.∆l[nDL], g3.∆l[nPR], isbloch, reorder=false)

    @test issymmetric(Mξ)

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mξ3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "param3d2mat"

@testset "symmetry of param3d2mat, nonuniform grid" begin
    # Create a grid.
    L₀ = 1e-9
    unit = PhysUnit(L₀)
    Npml = ([0,0,0], [0,0,0])
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5])
    g3 = Grid(unit, lprim, Npml, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = EncodedMaterial(PRIM, Material("Vacuum", ε=εvac))

    εdiel = 2.0
    diel = EncodedMaterial(PRIM, Material("Dielectric", ε=εdiel))

    # Create objects.
    dom_vac = Object(Box(g3.bounds), vac)
    obj_diel = Object(Box([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Sphere([0,0,0], 1), diel)

    # Add objects.
    ovec = Object3[]
    paramset = (SMat3Complex[], SMat3Complex[])
    add!(ovec, paramset, dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    param3d = create_param3d(N)
    obj3d = create_n3d(Object3, N)
    pind3d = create_n3d(ParamInd, N)
    oind3d = create_n3d(ObjInd, N)

    assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.isbloch)
    # Test the sanity the assigned param3d here.  It is relatively easy, and it was very helpful.

    smooth_param!(param3d, obj3d, pind3d, oind3d, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    Mξ = param3d2mat(param3d[nPR], PRIM, N, g3.∆l[nDL], g3.∆l[nPR], isbloch, reorder=false)

    # Ml
    ∆lprim = g3.∆l[nPR]
    ∆ldual = g3.∆l[nDL]

    Nx, Ny, Nz = N
    ∆X = repeat(∆ldual[nX], outer=(1,Ny,Nz))
    ∆Y = repeat(∆ldual[nY].', outer=(Nx,1,Nz))
    ∆Z = repeat(reshape(∆ldual[nZ], (1,1,Nz)), outer=(Nx,Ny,1))

    ∆L = [∆X[:].'; ∆Y[:].'; ∆Z[:].']
    ∆l = ∆L[:]
    Ml⁻¹ = spdiagm(1./∆l)

    # Ma
    ∆X = repeat(∆lprim[nX], outer=(1,Ny,Nz))
    ∆Y = repeat(∆lprim[nY].', outer=(Nx,1,Nz))
    ∆Z = repeat(reshape(∆lprim[nZ], (1,1,Nz)), outer=(Nx,Ny,1))

    ∆YZ = ∆Y .* ∆Z
    ∆ZX = ∆Z .* ∆X
    ∆XY = ∆X .* ∆Y

    ∆A = [∆YZ[:].'; ∆ZX[:].'; ∆XY[:].']
    ∆a = ∆A[:]
    Ma = spdiagm(∆a)

    Mξ_sym = Ma * Mξ * Ml⁻¹

    @test Mξ_sym ≈ Mξ_sym.'

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mξ3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "param3d2mat"

# Need to test if PML still works after the field-averaging operator applied to the material
# parameters, because effectively it makes the material inhomogeneous along the direction
# normal to the PML interface.

end  # @testset "param"
