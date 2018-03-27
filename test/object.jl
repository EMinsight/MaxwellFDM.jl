@testset "object" begin

const rtol = Base.rtoldefault(Float64)
const one⁻ = 1 - rtol  # slightly less than 1
const intv1⁻ = (-one⁻, one⁻)

@testset "Interval" begin
    vac = Material("vacuum")
    ge = PRIM
    evac = EncodedMaterial(ge, vac)

    box_vac = Object(Box([0,1], [2,3]), evac, [0.1, 0.15])
    b = bounds(box_vac)
    oi = OpenInterval(box_vac, nX)
    ci = ClosedInterval(box_vac, nY)

    @test bounds(oi) == (-1,1)
    @test bounds(ci) == (-0.5, 2.5)
    @test -1∉oi && 1∉oi && 0∈oi
    @test -0.5∈ci && 2.5∈ci && 1∈ci
    @test max∆l(oi) == 0.1
    @test max∆l(ci) == 0.15
    @test length(oi) == 2
    @test length(ci) == 3

    # ∆lmax = (bound[2]-bound[1]) / 10
    # center = mean(bound)
    # i1 = Interval1D(bound, ∆lmax)
    #
    # @test contains(i1, center)
    # @test all(contains.(i1, bound))
    # @test !contains(i1, bound[1]-eps())
    # @test center_(i1) ≈ center
end


# @testset "Interval1D" begin
#     bound = (sort(rand(2))...)
#     ∆lmax = (bound[2]-bound[1]) / 10
#     center = mean(bound)
#     i1 = Interval1D(bound, ∆lmax)
#
#     @test contains(i1, center)
#     @test all(contains.(i1, bound))
#     @test !contains(i1, bound[1]-eps())
#     @test center_(i1) ≈ center
# end
#
# @testset "Interval3D" begin
#     ba = sort(rand(3,2), 2)
#     bound = ([(ba[i,:]...) for i = 1:3]...)
#     ∆lmax = (/).((x->-(x...)).(bound), -10)
#     center = mean.(bound)
#     i3 = Interval3D(bound, ∆lmax)
#
#     @test contains(i3, center)
#     @test all(contains.(i3, [first.(bound), last.(bound)]))
#     @test center_(i3) ≈ center
# end

@testset "Object" begin
    mat = Material("material", ε=rand(), μ=rand())
    ge = PRIM
    emat = EncodedMaterial(ge, mat)

    @testset "Box" begin
        ba = sort(rand(3,2), 2)  # array
        b = (ba[:,1], ba[:,2])

        box_emat = Object(Box(b), emat)

        @test bounds(box_emat) ≈ b
        @test all(b .∈ box_emat)
        @test max∆l(box_emat) == fill(Inf,3)
        @test matparam(box_emat,PRIM)==emat.param[nPR] && matparam(box_emat,DUAL)==emat.param[nDL]
    end  # @testset "Box"

    @testset "Ellipsoid" begin
        c = rand(3)
        r = rand(3)
        b = (c-r, c+r)
        ∆lmax = r / 10

        el_emat = Object(Ellipsoid(c, r), emat, ∆lmax)

        @test bounds(el_emat) ≈ b
        @test all([(c′ = copy(c); c′[w] += s*r[w]; c′) for w = nXYZ, s = intv1⁻] .∈ el_emat)
        @test all([(c′ = copy(c); c′[nX] += sx*r[nX]; c′[nY] += sy*r[nY]; c′[nZ] += sz*r[nZ];
            all(bounds(el_emat)[nN] .≤ c′ .≤ bounds(el_emat)[nP])) for sx = intv1⁻, sy = intv1⁻, sz = intv1⁻])
        @test max∆l(el_emat) == ∆lmax
        @test matparam(el_emat,PRIM)==emat.param[nPR] && matparam(el_emat,DUAL)==emat.param[nDL]
    end  # @testset "Ellipsoid"

    @testset "Cylinder" begin
        c = rand(3)
        r = rand()
        h = rand()
        R = [r,r,h/2]
        a = [0,0,1]
        ∆lmax = R ./ 10

        cyl_emat = Object(Cylinder(c,r,a,h), emat, ∆lmax)

        b = (c-R, c+R)

        @test bounds(cyl_emat) ≈ b
        @test all([(c′ = copy(c); c′[w] += s*R[w]; c′) for w = (nX, nZ), s = intv1⁻] .∈ cyl_emat)
        @test all([(c′ = copy(c); c′[w] += s*R[w]; c′) for w = (nY, nZ), s = intv1⁻] .∈ cyl_emat)
        @test max∆l(cyl_emat) == ∆lmax
        @test matparam(cyl_emat,PRIM)==emat.param[nPR] && matparam(cyl_emat,DUAL)==emat.param[nDL]
    end  # @testset "Cylinder"

    @testset "Sphere" begin
        c = rand(3)
        r = rand()
        ∆lmax = r / 10

        sph_emat = Object(Sphere(c,r), emat, ∆lmax)

        R = [r,r,r]
        b = (c-R, c+R)

        @test bounds(sph_emat) ≈ b
        @test all([all([(c′ = copy(c); c′[w] += s*R[w]; c′) for s = intv1⁻] .∈ sph_emat) for w = nXYZ])
        @test max∆l(sph_emat) == fill(∆lmax,3)
        @test matparam(sph_emat,PRIM)==emat.param[nPR] && matparam(sph_emat,DUAL)==emat.param[nDL]
    end  # @testset "Sphere"
end  # @testset "Object"

@testset "Object Vector" begin
    ovec = Object3[]
    paramset = (SMat3Complex[], SMat3Complex[])
    vac = Material("vacuum")
    Si = Material("Si", ε = 12)

    ge = PRIM
    evac = EncodedMaterial(ge, vac)
    eSi = EncodedMaterial(ge, Si)

    box_evac = Object(Box((rand(3),rand(3))), evac)
    add!(ovec, paramset, box_evac)

    el_eSi = Object(Ellipsoid(rand(3),rand(3)), eSi)
    add!(ovec, paramset, el_eSi)

    box2_evac = Object(Box((rand(3),rand(3))), evac)
    add!(ovec, paramset, box2_evac)

    @test paramind(box_evac,PRIM)==1 && paramind(box_evac,DUAL)==1
    @test objind(box_evac) == 1

    @test paramind(el_eSi,PRIM)==2 && paramind(el_eSi,DUAL)==1
    @test objind(el_eSi) == 2

    @test paramind(box2_evac,PRIM)==1 && paramind(box2_evac,DUAL)==1
    @test objind(box2_evac) == 3
end  # @testset "Object vector"

@testset "periodize" begin
    # Square lattice
    Si = Material("Si", ε = 12)
    ge = PRIM
    eSi = EncodedMaterial(ge, Si)
    c_eSi = Object(Cylinder([0,0,0], 1, [0,0,1], 5), eSi)

    b = Box([0,0,0], [10,10,5])
    A = [1 0 0; 0 1 0; 0 0 5]'
    obj_array = periodize(c_eSi, A, b)

    @test length(obj_array) == (11-2)^2
    @test all(map(x->(x==max∆l(c_eSi)), max∆l.(obj_array)))
    @test all(map(x->(x==matparam(c_eSi,PRIM)), matparam.(obj_array,PRIM)))
    @test all(map(x->(x==matparam(c_eSi,DUAL)), matparam.(obj_array,DUAL)))
end  # @testset "periodize"


# @testset "union" begin
#     b1 = Box(((0,1),(0,1),(0,1)))
#     b2 = Box(((-1,0),(-1,0),(-1,0)))
#     s = sphere((0,0,0), 1)
#     lsf_union = union(lsf(b1), lsf(s), lsf(b2))
#     @test all(contains.(lsf_union, [(1,1,1), (-1,-1,-1), (-1/√3,-1/√3,1/√3), (1/√3,1/√3,-1/√3)]))
# end
#
# @testset "intersect" begin
#     b1 = Box(((-1,2),(-1,2),(-1,2)))
#     b2 = Box(((-2,1),(-2,1),(-2,1)))
#     lsf_intersect = intersect(lsf(b1), lsf(b2))
#     @test all(contains.(lsf_intersect, [(x,y,z) for x = -1:1, y = -1:1, z = -1:1]))
#     @test all(!contains.(lsf_intersect, [(x,y,z) for x = (-2,2), y = (-2,2), z = (-2,2)]))
# end
#
# @testset "flip" begin
#     b = Box(((0,1),(0,1),(0,1)))
#     lsf0 = lsf(b)
#     lsf_x = flip(lsf0, XX, 0)
#     lsf_y = flip(lsf0, YY, 0)
#     lsf_z = flip(lsf0, ZZ, 0)
#
#     @test all(contains.(lsf0, [(x,y,z) for x = 0:1, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_x, [(x,y,z) for x = -1:0, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_y, [(x,y,z) for x = 0:1, y = -1:0, z = 0:1]))
#     @test all(contains.(lsf_z, [(x,y,z) for x = 0:1, y = 0:1, z = -1:0]))
# end
#
# @testset "shift" begin
#     b = Box(((0,1),(0,1),(0,1)))
#     lsf0 = lsf(b)
#     lsf_x = shift(lsf0, XX, -1)
#     lsf_y = shift(lsf0, YY, -1)
#     lsf_z = shift(lsf0, ZZ, -1)
#
#     @test all(contains.(lsf0, [(x,y,z) for x = 0:1, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_x, [(x,y,z) for x = -1:0, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_y, [(x,y,z) for x = 0:1, y = -1:0, z = 0:1]))
#     @test all(contains.(lsf_z, [(x,y,z) for x = 0:1, y = 0:1, z = -1:0]))
# end

end  # @testset "object"
