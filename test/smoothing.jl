@testset "smoothing" begin

@testset "sort8!" begin
    v = rand(8)
    ind = collect(1:8)
    @inferred(MaxwellFD3D.sort8!(ind,v))
    @test issorted(v[ind])
    @test @inferred(MaxwellFD3D.countdiff(ind,v)) == (8,7)

    @inferred(MaxwellFD3D.sort8!(v))
    @test issorted(v)
    @test @inferred(MaxwellFD3D.countdiff(v)) == 8
end  # @testset "sort8!"

# @testset "Object" begin
#     vac = Material("Vacuum")
#     ivac = EncodedMaterial(PRIM, vac)
#     box = Box(((0,1), (0,1), (0,1)))
#     obj = Object(ivac, box)
#     obj_array = Object.(ivac, [box,box,box])  # vectorization over shapes
#     @test obj_array == [obj, obj, obj]
# end  # @testset "Object"
#
@testset "smoothing" begin
    # Create a grid.
    L₀ = 1e-9
    unit = PhysUnit(L₀)
    Npml = ([0,0,0], [0,0,0])
    N = [3,3,3]
    M = sum(Npml) + N
    ∆ldual = map(m->ones(m), M)  # tuple of arrays
    L = sum(∆ldual)
    l₀ = L / 2
    # lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-2.0, -1.0, 0.0, 1.0, 2.0])
    # ebc = [BLOCH, PPC, PDC]
    lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5])
    ebc = [BLOCH, BLOCH, BLOCH]
    g3 = Grid(unit, lprim, Npml, ebc)

    # Create materials.
    εvac = 1.0
    vac = EncodedMaterial(PRIM, Material("Vacuum", ε=εvac))

    εdiel = 2.0
    diel = EncodedMaterial(PRIM, Material("Dielectric", ε=εdiel))

    # Create objects.
    dom_vac = setmat!(Box(g3.bounds), vac)
    obj_diel = setmat!(Box([0,0,0], [1,1,1]), diel)
    # obj_diel = setmat!(Sphere([0,0,0], 1), diel)

    # Add objects.
    ovec = Object3[]
    paramset = (CMatrix3[], CMatrix3[])
    add!(ovec, paramset, dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    param3d = create_param3d_zeros(g3.N)
    obj3d = create_n3d_zeros(Object3, g3.N)
    pind3d = create_n3d_zeros(ParamInd, g3.N)
    oind3d = create_n3d_zeros(ObjInd, g3.N)

    assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.ebc)
    smooth_param!(param3d, obj3d, pind3d, oind3d, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    ε3d = view(param3d[nPR], 1:M[nX], 1:M[nY], 1:M[nZ], 1:3, 1:3)

    # Construct an expected ε3d.
    ε3dexp = Array{Complex128}(3,3,3,3,3)
    rvol = 0.5  # all nonzero rvol used in this test is 0.5
    εh = 1 / (rvol/εdiel + (1-rvol)/εvac)  # harmonic average
    εa = rvol*εdiel + (1-rvol)*εvac  # arithmetic average

    # Initialize ε3dexp.
    I = eye(3)
    for k = 1:3, j = 1:3, i = 1:3
        ε3dexp[i,j,k,:,:] = εvac * I
    end

    # Yee's cell at (2,2,2)
    nout = normalize([-1,-1,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,2,2,:,:] = εsm
    nout = normalize([0,-1,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,2,2,1,1] = εsm[1,1]
    nout = normalize([-1,0,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,2,2,2,2] = εsm[2,2]
    nout = normalize([-1,-1,0]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,2,2,3,3] = εsm[3,3]

    # Yee's cell at (3,2,2)
    nout = normalize([1,-1,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,2,2,:,:] = εsm
    ε3dexp[3,2,2,1,1] = εvac
    nout = normalize([1,0,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,2,2,2,2] = εsm[2,2]
    nout = normalize([1,-1,0]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,2,2,3,3] = εsm[3,3]

    # Yee's cell at (2,3,2)
    nout = normalize([-1,1,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,3,2,:,:] = εsm
    nout = normalize([0,1,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,3,2,1,1] = εsm[1,1]
    ε3dexp[2,3,2,2,2] = εvac
    nout = normalize([-1,1,0]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,3,2,3,3] = εsm[3,3]

    # Yee's cell at (2,2,3)
    nout = normalize([-1,-1,1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,2,3,:,:] = εsm
    nout = normalize([0,-1,1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,2,3,1,1] = εsm[1,1]
    nout = normalize([-1,0,1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,2,3,2,2] = εsm[2,2]
    ε3dexp[2,2,3,3,3] = εvac

    # Yee's cell at (2,3,3)
    nout = normalize([-1,1,1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,3,3,:,:] = εsm
    nout = normalize([0,1,1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[2,3,3,1,1] = εsm[1,1]
    ε3dexp[2,3,3,2,2] = εvac
    ε3dexp[2,3,3,3,3] = εvac

    # Yee's cell at (3,2,3)
    nout = normalize([1,-1,1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,2,3,:,:] = εsm
    ε3dexp[3,2,3,1,1] = εvac
    nout = normalize([-1,0,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,2,3,2,2] = εsm[2,2]
    ε3dexp[3,2,3,3,3] = εvac

    # Yee's cell at (3,3,2)
    nout = normalize([1,1,-1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,3,2,:,:] = εsm
    ε3dexp[3,3,2,1,1] = εvac
    ε3dexp[3,3,2,2,2] = εvac
    nout = normalize([1,1,0]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,3,2,3,3] = εsm[3,3]

    # Yee's cell at (3,3,3)
    nout = normalize([1,1,1]); P = nout*nout'; εsm = εh * P + εa * (I-P); ε3dexp[3,3,3,:,:] = εsm
    ε3dexp[3,3,3,1,1] = εvac
    ε3dexp[3,3,3,2,2] = εvac
    ε3dexp[3,3,3,3,3] = εvac

    @test ε3d ≈ ε3dexp
end

end  # @testset "object"
