using MaxwellFDM
using StaticArrays
using GeometryPrimitives
using Base.Test

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellFDM" begin

include("enumtype.jl")
include("base.jl")
include("phys.jl")
include("grid.jl")
include("object.jl")
include("gridgen.jl")
include("material.jl")
include("smoothing.jl")
include("operator.jl")
include("pml.jl")
# include("source.jl")

# end  # @testset "MaxwellFDM"
