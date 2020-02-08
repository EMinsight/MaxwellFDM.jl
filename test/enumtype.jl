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

@testset "ft2gt" begin
    @test MaxwellFDM.ft2gt(EE,EE) == PRIM
    @test MaxwellFDM.ft2gt(HH,EE) == DUAL
    @test MaxwellFDM.ft2gt(EE,HH) == DUAL
    @test MaxwellFDM.ft2gt(HH,HH) == PRIM
end

@testset "gt_Fw" begin
    # 3D
    @test MaxwellFDM.gt_Fw(1, EE, SVector(EE,EE,EE)) == SVector(DUAL,PRIM,PRIM)
    @test MaxwellFDM.gt_Fw(2, EE, SVector(EE,EE,EE)) == SVector(PRIM,DUAL,PRIM)
    @test MaxwellFDM.gt_Fw(3, EE, SVector(EE,EE,EE)) == SVector(PRIM,PRIM,DUAL)
    @test MaxwellFDM.gt_Fw(4, EE, SVector(EE,EE,EE)) == SVector(PRIM,PRIM,PRIM)

    # 2D
    @test MaxwellFDM.gt_Fw(1, EE, SVector(EE,EE)) == SVector(DUAL,PRIM)
    @test MaxwellFDM.gt_Fw(2, EE, SVector(EE,EE)) == SVector(PRIM,DUAL)
    @test MaxwellFDM.gt_Fw(3, EE, SVector(EE,EE)) == SVector(PRIM,PRIM)
end

end  # @testset "enumtype"
