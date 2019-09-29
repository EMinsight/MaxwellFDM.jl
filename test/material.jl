@testset "material" begin

@testset "Material" begin
    m_nn = Material("(num,num)", ε = 2, μ = 3)
    m_vv = Material("(vec,vec)", ε = [2.1,2.2,2.3], μ = [3.1,3.2,3.3])
    m_mm = Material("(mat,mat)", ε = [2 2 2; 2 2 2; 2 2 2], μ = [3 3 3; 3 3 3; 3 3 3])

    @test matparam(m_nn,EE) == diagm(0=>[2,2,2]) && matparam(m_nn,HH) == diagm(0=>[3,3,3]) && string(m_nn) == "(num,num)"
    @test matparam(m_vv,EE) == diagm(0=>[2.1,2.2,2.3]) && matparam(m_vv,HH) == diagm(0=>[3.1,3.2,3.3]) && string(m_vv) == "(vec,vec)"
    @test matparam(m_mm,EE) == fill(2, (3,3)) && matparam(m_mm,HH) == fill(3, (3,3)) && string(m_mm) == "(mat,mat)"
end

end  # @testset "material"
