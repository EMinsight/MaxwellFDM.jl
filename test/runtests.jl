using MaxwellWave
using Test
using Statistics: mean
using StaticArrays

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellWave" begin

include("source.jl")

# end  # @testset "MaxwellWave"
