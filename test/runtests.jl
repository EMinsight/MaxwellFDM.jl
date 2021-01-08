using MaxwellWave
using Test
using Statistics: mean
using StaticArrays

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellWave" begin

include("enumtype.jl")
include("util.jl")
include("phys.jl")
include("material.jl")
include("object.jl")
include("smoothing.jl")
include("param.jl")
include("source.jl")

# end  # @testset "MaxwellWave"
