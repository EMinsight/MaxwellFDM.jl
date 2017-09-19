@testset "enumtype" begin

@testset "instances" begin
    @test instances(Axis) == (XX, YY, ZZ)
    @test instances(Dir) == (HRZ, VRT)
    @test instances(GridType) == (PRIM, DUAL)
    @test instances(FieldType) == (EE, HH)
    @test instances(Sign) == (NEG, POS)
    @test instances(BC) == (PERIODIC, PEC, PMC)
    @test instances(EBC) == (BLOCH, PPC, PDC)
    @test instances(PML) == (SCPML, UPML)
end

@testset "numel" begin
    @test numel(Axis) == 3
    @test numel(Dir) == 2
    @test numel(GridType) == 2
    @test numel(FieldType) == 2
    @test numel(Sign) == 2
    @test numel(BC) == 3
    @test numel(EBC) == 3
    @test numel(PML) == 2
end

@testset "integers" begin
    @test Int.(SVector(instances(Axis))) == nXYZ == [nX, nY, nZ]
    @test Int.(SVector(instances(Dir))) == nHV == [nHRZ, nVRT]
    @test Int.(SVector(instances(Sign))) == nNP == [nN, nP]
    @test Int.(SVector(instances(GridType))) == nPD == [nPR, nDL]
end

@testset "next and alter" begin
    @test next3(XX) == [YY, ZZ, XX]
    @test next3(YY) == [ZZ, XX, YY]
    @test next3(ZZ) == [XX, YY, ZZ]

    @test next2(XX) == [YY, ZZ]
    @test next2(YY) == [ZZ, XX]
    @test next2(ZZ) == [XX, YY]

    @test next1(XX) == YY
    @test next1(YY) == ZZ
    @test next1(ZZ) == XX

    @test prev3(XX) == [ZZ, YY, XX]
    @test prev3(YY) == [XX, ZZ, YY]
    @test prev3(ZZ) == [YY, XX, ZZ]

    @test prev2(XX) == [ZZ, YY]
    @test prev2(YY) == [XX, ZZ]
    @test prev2(ZZ) == [YY, XX]

    @test prev1(XX) == ZZ
    @test prev1(YY) == XX
    @test prev1(ZZ) == YY

    @test alter(HRZ) == VRT
    @test alter(VRT) == HRZ

    @test alter(PRIM) == DUAL
    @test alter(DUAL) == PRIM

    @test alter(EE) == HH
    @test alter(HH) == EE

    @test alter(NEG) == POS
    @test alter(POS) == NEG

    @test alter(SCPML) == UPML
    @test alter(UPML) == SCPML
end

end  # @testset "enumtype"
