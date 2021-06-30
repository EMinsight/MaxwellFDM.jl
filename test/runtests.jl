using MaxwellFDFD
using Test
using Statistics: mean
using AbbreviatedTypes

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellFDFD" begin

include("source.jl")

# end  # @testset "MaxwellFDFD"
