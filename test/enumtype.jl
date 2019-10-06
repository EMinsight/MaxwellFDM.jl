@testset "enumtype" begin

@testset "instances" begin
    @test instances(FieldType) == (EE, HH)
    @test instances(BC) == (PERIODIC, CONDUCTING)
end

@testset "numel" begin
    @test numel(FieldType) == 2
    @test numel(BC) == 2
end

@testset "integers" begin
    @test Int.(SVector(instances(FieldType))) == nEH == [nE, nH]
end

@testset "alter" begin
    @test alter(EE) == HH
    @test alter(HH) == EE

    @test alter(nE) == nH
    @test alter(nH) == nE
end

end  # @testset "enumtype"
