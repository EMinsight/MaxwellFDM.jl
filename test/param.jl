@testset "param" begin

# Need to test symmetry of param3dmat for different boundary conditions.
@testset "symmetry of param_arr2mat, uniform grid" begin
    # Create a grid.
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5])
    g3 = Grid(lprim, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = Material("Vacuum", ε=εvac)

    εdiel = 2.0
    diel = Material("Dielectric", ε=εdiel)

    # Create objects.
    dom_vac = Object(Box(g3.bounds), vac)
    obj_diel = Object(Box([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Sphere([0,0,0], 1), diel)

    # Add objects.
    ovec = Object{3}[]
    paramset = (SSComplex3[], SSComplex3[])
    add!(ovec, paramset, dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    ε3d = create_param_array(N)
    εobj3d = create_p_storage(Object{3}, N)
    εind3d = create_p_storage(ParamInd, N)
    εoind3d = create_p_storage(ObjInd, N)

    μ3d = create_param_array(N)
    μobj3d = create_p_storage(Object{3}, N)
    μind3d = create_p_storage(ParamInd, N)
    μoind3d = create_p_storage(ObjInd, N)

    boundft = SVector(EE,EE,EE)
    assign_param!((ε3d,μ3d), (εobj3d,μobj3d), (εind3d,μind3d), (εoind3d,μoind3d), boundft, ovec, g3.ghosted.τl, g3.isbloch)
    # Test the sanity the assigned param3d here.  It is relatively easy, and it was very helpful.

    ft = EE
    smooth_param!(ε3d, εobj3d, εind3d, εoind3d, ft, boundft, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    Mε = param_arr2mat(ε3d, ft, boundft, N, g3.∆l[nDL], g3.∆l[nPR], g3.isbloch, reorder=false)

    @test issymmetric(Mε)

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mε3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "param_arr2mat"

@testset "symmetry of param_arr2mat, nonuniform grid" begin
    # Create a grid.
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5])
    g3 = Grid(lprim, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = Material("Vacuum", ε=εvac)

    εdiel = 2.0
    diel = Material("Dielectric", ε=εdiel)

    # Create objects.
    dom_vac = Object(Box(g3.bounds), vac)
    obj_diel = Object(Box([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Sphere([0,0,0], 1), diel)

    # Add objects.
    ovec = Object{3}[]
    paramset = (SSComplex3[], SSComplex3[])
    add!(ovec, paramset, dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    ε3d = create_param_array(N)
    εobj3d = create_p_storage(Object{3}, N)
    εind3d = create_p_storage(ParamInd, N)
    εoind3d = create_p_storage(ObjInd, N)

    μ3d = create_param_array(N)
    μobj3d = create_p_storage(Object{3}, N)
    μind3d = create_p_storage(ParamInd, N)
    μoind3d = create_p_storage(ObjInd, N)

    boundft = SVector(EE,EE,EE)
    assign_param!((ε3d,μ3d), (εobj3d,μobj3d), (εind3d,μind3d), (εoind3d,μoind3d), boundft, ovec, g3.ghosted.τl, g3.isbloch)
    # Test the sanity the assigned param3d here.  It is relatively easy, and it was very helpful.

    ft = EE
    smooth_param!(ε3d, εobj3d, εind3d, εoind3d, ft, boundft, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    Mε = param_arr2mat(ε3d, ft, boundft, N, g3.∆l[nDL], g3.∆l[nPR], g3.isbloch, reorder=false)

    # Ml
    ∆lprim = g3.∆l[nPR]
    ∆ldual = g3.∆l[nDL]

    Nx, Ny, Nz = N
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
end  # @testset "param_arr2mat"

# Need to test if PML still works after the field-averaging operator applied to the material
# parameters, because effectively it makes the material inhomogeneous along the direction
# normal to the PML interface.

end  # @testset "param"
