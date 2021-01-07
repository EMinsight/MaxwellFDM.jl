@testset "object" begin

@testset "Object" begin
    mat = Material{3,3}("material", ε=rand(), μ=rand())

    @testset "Box" begin
        ba = sort(rand(3,2), dims=2)  # array
        b = (ba[:,1], ba[:,2])

        box_mat = Object(Box(b), mat)  # new convenience constructor

        @test shape(box_mat) == Box(b)
        @test max∆l(box_mat) == fill(Inf,3)
        @test matparam(box_mat,EE)==mat.param[nE] && matparam(box_mat,HH)==mat.param[nH]
    end  # @testset "Box"
end  # @testset "Object"

# @testset "periodize" begin
#     # Square lattice
#     Si = Material("Si", ε = 12)
#     ge = PRIM
#     eSi = EncodedMaterial(ge, Si)
#     c_eSi = Object(Cylinder([0,0,0], 1, 5, [0,0,1]), eSi)
#
#     b = Box([0,0,0], [10,10,5])
#     A = [1 0 0; 0 1 0; 0 0 5]'
#     obj_array = periodize(c_eSi, A, b)
#
#     @test length(obj_array) == 11^2
#     @test all(map(x->(x==max∆l(c_eSi)), max∆l.(obj_array)))
#     @test all(map(x->(x==matparam(c_eSi,PRIM)), matparam.(obj_array,PRIM)))
#     @test all(map(x->(x==matparam(c_eSi,DUAL)), matparam.(obj_array,DUAL)))
#
#     oind2obj = Object{3}[]
#     pind2matprm = (SSComplex3[], SSComplex3[])
#     add!(oind2obj, pind2matprm, obj_array)
# end  # @testset "periodize"


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
#     lsf_x = flip(lsf0, X̂, 0)
#     lsf_y = flip(lsf0, Ŷ, 0)
#     lsf_z = flip(lsf0, Ẑ, 0)
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
#     lsf_x = shift(lsf0, X̂, -1)
#     lsf_y = shift(lsf0, Ŷ, -1)
#     lsf_z = shift(lsf0, Ẑ, -1)
#
#     @test all(contains.(lsf0, [(x,y,z) for x = 0:1, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_x, [(x,y,z) for x = -1:0, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_y, [(x,y,z) for x = 0:1, y = -1:0, z = 0:1]))
#     @test all(contains.(lsf_z, [(x,y,z) for x = 0:1, y = 0:1, z = -1:0]))
# end

end  # @testset "object"
