@testset "material" begin

@testset "Material" begin
    m_nn = Material("(num,num)", ε = 2, μ = 3)
    m_vv = Material("(vec,vec)", ε = [2.1,2.2,2.3], μ = [3.1,3.2,3.3])
    m_mm = Material("(mat,mat)", ε = [2 2 2; 2 2 2; 2 2 2], μ = [3 3 3; 3 3 3; 3 3 3])

    @test m_nn.ε == diagm([2,2,2]) && m_nn.μ == diagm([3,3,3]) && string(m_nn) == "(num,num)"
    @test m_vv.ε == diagm([2.1,2.2,2.3]) && m_vv.μ == diagm([3.1,3.2,3.3]) && string(m_vv) == "(vec,vec)"
    @test m_mm.ε == fill(2, (3,3)) && m_mm.μ == fill(3, (3,3)) && string(m_mm) == "(mat,mat)"
end  # @testset "material"

@testset "EncodedMaterial" begin
    si = Material("Si", ε = 12)
    eSi = EncodedMaterial(DUAL, si)

    @test matparam(eSi, PRIM) == eye(3)
    @test matparam(eSi, DUAL) == diagm([12,12,12])
    @test string(eSi) == "Si"
end

end  # @testset "material"
