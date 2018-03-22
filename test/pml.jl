@testset "pml" begin

@testset "PMLParam" begin
    pml = PMLParam()
    @test pml.m===4.0 && pml.R===exp(-16) && pml.κmax===1.0 && pml.amax===0.0 && pml.ma===4.0
end  # @testset "PMLParam"

end  # @testset "pml"
