@testset "util" begin

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
    @test @inferred(MaxwellFDM.t_ind(t23, CartesianIndex(1,2,1))) == (1,5,3)
    @test @inferred(MaxwellFDM.t_ind(t23, (PRIM,DUAL,PRIM))) == (1,5,3)
    @test @inferred(MaxwellFDM.t_ind(t23, SVector(PRIM,DUAL,PRIM))) == (1,5,3)

    @test @inferred(MaxwellFDM.t_ind(t2a, 2, 3)) == [2,6]
    @test @inferred(MaxwellFDM.t_ind(t2sa, NEG, POS, NEG)) == [1,5,3]
    @test @inferred(MaxwellFDM.t_ind(t2sa, (NEG,POS,NEG))) == [1,5,3]
    @test @inferred(MaxwellFDM.t_ind(t2sa, SVector(NEG,POS,NEG))) == [1,5,3]

    @test @inferred(MaxwellFDM.t_ind(t3a, 2, 3, 4)) == [2,7,12]
    @test @inferred(MaxwellFDM.t_ind(t3a, CartesianIndex(2,3,4))) == [2,7,12]
    @test @inferred(MaxwellFDM.t_ind(t3a, SVector(2,3,4))) == [2,7,12]
end

@testset "sort_ind!" begin
    v = MArray{Tuple{2,2,2}}(rand(8))
    ind = MArray{Tuple{2,2,2}}(collect(1:8))
    @inferred(MaxwellFDM.sort_ind!(ind,v))
    @test issorted(v[ind])
    @test @inferred(MaxwellFDM.countdiff(ind,v)) == (8,8)
end  # @testset "sort8!"

end  # @testset "util"
