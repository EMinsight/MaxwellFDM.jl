using MaxwellFDM
using Test
using Statistics: mean
using LinearAlgebra, SparseArrays, StaticArrays, StaggeredGridCalculus, GeometryPrimitives

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellFDM" begin

include("enumtype.jl")
include("util.jl")
include("phys.jl")
include("object.jl")
include("material.jl")
include("smoothing.jl")
include("param.jl")
include("source.jl")

# end  # @testset "MaxwellFDM"
