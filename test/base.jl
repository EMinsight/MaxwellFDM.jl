@testset "base" begin

@testset "getindex" begin
    t3 = (1,2,3)
    t32 = ((1,2), (3,4), (5,6))
    t23 = ((1,2,3), (4,5,6))
    t2sa = (SVector(1,2,3), SVector(4,5,6))
    t2a = ([1,2,3], [4,5,6])
    t3a = ([1,2,3,4], [5,6,7,8], [9,10,11,12])

    # @test t_ind(t3, (2,3,1)) == (2,3,1)
    # @test t_ind(t32, 1, 2, 1) == (1,4,5)
    @test @inferred(MaxwellFDM.t_ind(t23, PRIM, DUAL, PRIM)) == (1,5,3)
    @test @inferred(MaxwellFDM.t_ind(t23, (PRIM,DUAL,PRIM))) == (1,5,3)
    @test @inferred(MaxwellFDM.t_ind(t23, SVector(PRIM,DUAL,PRIM))) == (1,5,3)
    @test @inferred(MaxwellFDM.t_ind(t2a, 2, 3)) == [2,6]
    @test @inferred(MaxwellFDM.t_ind(t2sa, NEG, POS, NEG)) == [1,5,3]
    @test @inferred(MaxwellFDM.t_ind(t2sa, (NEG,POS,NEG))) == [1,5,3]
    @test @inferred(MaxwellFDM.t_ind(t2sa, SVector(NEG,POS,NEG))) == [1,5,3]
    @test @inferred(MaxwellFDM.t_ind(t3a, 2, 3, 4)) == [2,7,12]
    @test @inferred(MaxwellFDM.t_ind(t3a, SVector(2,3,4))) == [2,7,12]
    @test @inferred(MaxwellFDM.t_ind(t3a, ((1,2,3), (2,3,4)))) == ([1,6,11], [2,7,12])
    @test @inferred(MaxwellFDM.t_ind(t3a, (SVector(1,2,3), SVector(2,3,4)))) == ([1,6,11], [2,7,12])
end

# @testset "dot" begin
#     t1 = (1,2,3)
#     t2 = (4,5,6)
#     @test t1⋅t2 == sum(t1 .* t2)
# end

# @testset "isapprox" begin
#     aa = [[1,2], [3,4]]
#     tt = ((1,2), (3,4))
#     at = [(1,2), (3,4)]
#     ta = ([1,2], [3,4])
#
#     @test aa ≈ aa
#     @test tt ≈ tt
#     @test_broken at ≈ at
#     @test ta ≈ ta
#     @test aa ≉ tt
#     @test tt ≉ at
#     @test at ≉ ta
#     @test_throws MethodError aa ≉ at
#     @test tt ≉ ta
#     @test aa ≉ ta
# end

@testset "newtsol" begin
    f(x) = x^2-1
    f′(x) = 2x
    @test ((xsol,isconverged) = MaxwellFDM.newtsol(2., f, f′); isconverged && xsol ≈ 1)
    @test ((xsol,isconverged) = MaxwellFDM.newtsol(2., f); isconverged && xsol ≈ 1)
end

end  # @testset "base"
