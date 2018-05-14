@testset "operator" begin

@testset "create_∂" begin
    N = SVector(8,9,10)
    # N = SVector(3,3,3)
    M = prod(N)
    for nw = nXYZ
    # for nw = (1,)
        Nw = N[nw]
        sub′ = Vector{Int}(3)

        for ns = (-1,1)
        # for ns = (1,)
            ∂ws = spzeros(M,M)

            for ind = 1:M
                ∂ws[ind,ind] = -ns  # diagonal entries

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
                ∂ws[ind, ind′] += ns  # off-diagonal entry
            end
            @test create_∂(nw, ns, N) == ∂ws
        end
    end
end  # @testset "create_∂"

N = SVector(3,4,5)
M = prod(N)
r = reshape(collect(1:3M), M, 3)'[:]  # index mapping from block matrix to narrowly banded matrix
Z = spzeros(M,M)

@testset "create_curl for U" begin
    # Construct Cu for a uniform grid and BLOCH boundaries.
    Cu = create_curl(PRIM, N, reorder=false)

    # Examine the overall coefficients.
    @test all(any(Cu.≠0, 1))  # no zero columns
    @test all(any(Cu.≠0, 2))  # no zero rows
    @test all(sum(Cu, 2) .== 0)  # all row sums are zero, because Cu * ones(M) = 0

    ∂x = (nw = 1; create_∂(nw, 1, N))
    ∂y = (nw = 2; create_∂(nw, 1, N))
    ∂z = (nw = 3; create_∂(nw, 1, N))
    @test Cu == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Construct Cu for a nonuniform grid and general boundaries.
    ∆ldual = rand.(N.data)
    ebc = SVector(BLOCH, PPC, PDC)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    Cu = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine Cu.
    ∂x = (nw = 1; create_∂(nw, 1, N, ∆ldual[nw], ebc[nw], e⁻ⁱᵏᴸ[nw]))
    ∂y = (nw = 2; create_∂(nw, 1, N, ∆ldual[nw], ebc[nw], e⁻ⁱᵏᴸ[nw]))
    ∂z = (nw = 3; create_∂(nw, 1, N, ∆ldual[nw], ebc[nw], e⁻ⁱᵏᴸ[nw]))
    @test Cu == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Examine reordering.
    Cu_reorder = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=true)
    @test Cu_reorder == Cu[r,r]
end  # @testset "create_curl for U"

@testset "create_curl for V" begin
    # Construct Cv for a uniform grid and BLOCH boundaries.
    Cv = create_curl(DUAL, N, reorder=false)

    # Examine the overall coefficients.
    @test all(any(Cv.≠0, 1))  # no zero columns
    @test all(any(Cv.≠0, 2))  # no zero rows
    @test all(sum(Cv, 2) .== 0)  # all row sums are zero, because Cv * ones(sum(Min)) = 0

    ∂x = (nw = 1; create_∂(nw, -1, N))
    ∂y = (nw = 2; create_∂(nw, -1, N))
    ∂z = (nw = 3; create_∂(nw, -1, N))
    @test Cv == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Construct Cv for a nonuniform grid and general boundaries.
    ∆lprim = rand.(N.data)
    ebc = SVector(BLOCH, PPC, PDC)
    e⁻ⁱᵏᴸ = @SVector rand(Complex128, 3)

    Cv = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Examine Cv.
    ∂x = (nw = 1; create_∂(nw, -1, N, ∆lprim[nw], ebc[nw], e⁻ⁱᵏᴸ[nw]))
    ∂y = (nw = 2; create_∂(nw, -1, N, ∆lprim[nw], ebc[nw], e⁻ⁱᵏᴸ[nw]))
    ∂z = (nw = 3; create_∂(nw, -1, N, ∆lprim[nw], ebc[nw], e⁻ⁱᵏᴸ[nw]))
    @test Cv == [Z -∂z ∂y;
                 ∂z Z -∂x;
                 -∂y ∂x Z]

    # Examine reordering
    Cv_reorder = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=true)
    @test Cv_reorder == Cv[r,r]
end  # @testset "create_curl for V"

@testset "curl of curl" begin
    # Construct Cu and Cv for a uniform grid and BLOCH boundaries.
    ∆ldual = ones.(N.data)
    ∆lprim = ones.(N.data)
    ebc =  SVector(BLOCH, PPC, PDC)
    e⁻ⁱᵏᴸ = @SVector ones(3)

    Cu = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)
    Cv = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)

    # Test symmetry of each block.
    for i = nXYZ
        for j = next2(i)
            -Cv[(i-1)*M+1:i*M,(j-1)*M+1:j*M]' == Cu[(i-1)*M+1:i*M,(j-1)*M+1:j*M]
        end
    end

    # Construct Cv * Cu for all BLOCH.
    ebc =  @SVector fill(BLOCH, 3)
    Cu = create_curl(PRIM, N, ∆ldual, ebc, e⁻ⁱᵏᴸ, reorder=false)
    Cv = create_curl(DUAL, N, ∆lprim, ebc, e⁻ⁱᵏᴸ, reorder=false)
    A = Cv * Cu

    # Test curl of curl.
    @test all(diag(A) .== 4)  # all diagonal entries are 4
    @test all(sum(A.≠0, 2) .== 13)  # 13 nonzero entries per row
    @test A == A'  # Hermitian

    B = A - 4*speye(A)
    @test all(abs.(B[B.≠0]).==1)  # all nonzero off-diagonal entries are ±1
end  # @testset "curl of curl"

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

# Need to test symmetry of param3dmat for different boundary conditions.
@testset "symmetry of param3d2mat, uniform grid" begin
    using MaxwellFDM, GeometryPrimitives, StaticArrays

    # Create a grid.
    L₀ = 1e-9
    unit = PhysUnit(L₀)
    Npml = ([0,0,0], [0,0,0])
    ebc = [BLOCH, BLOCH, BLOCH]
    lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5])
    # ebc = [BLOCH, PPC, PDC]
    # lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-2.0, -1.0, 0.0, 1.0, 2.0])
    g3 = Grid(unit, lprim, Npml, ebc)
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

    assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.ebc)
    # Test the sanity the assigned param3d here.  It is relatively easy, and it was very helpful.

    smooth_param!(param3d, obj3d, pind3d, oind3d, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    Mξ = param3d2mat(param3d[nPR], PRIM, N, g3.∆l[nDL], g3.∆l[nPR], ebc, reorder=false)

    @test issymmetric(Mξ)

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mξ3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "param3d2mat"

@testset "symmetry of param3d2mat, nonuniform grid" begin
    using MaxwellFDM, GeometryPrimitives, StaticArrays

    # Create a grid.
    L₀ = 1e-9
    unit = PhysUnit(L₀)
    Npml = ([0,0,0], [0,0,0])
    ebc = [BLOCH, BLOCH, BLOCH]
    lprim = ([-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5])
    # ebc = [BLOCH, PPC, PDC]
    # lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-2.0, -1.0, 0.0, 1.0, 2.0])
    g3 = Grid(unit, lprim, Npml, ebc)
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

    assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.ebc)
    # Test the sanity the assigned param3d here.  It is relatively easy, and it was very helpful.

    smooth_param!(param3d, obj3d, pind3d, oind3d, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

    Mξ = param3d2mat(param3d[nPR], PRIM, N, g3.∆l[nDL], g3.∆l[nPR], ebc, reorder=false)

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

end  # @testset "operator"
