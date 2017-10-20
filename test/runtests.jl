using MaxwellFD3D
using StaticArrays
using Base.Test

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellFD3D" begin

# include("enumtype.jl")
# include("base.jl")
# include("phys.jl")
# include("grid.jl")
# include("shape.jl")
# include("gridgen.jl")
# include("material.jl")
# include("smoothing.jl")
include("operator.jl")
# include("pml.jl")

# end  # @testset "MaxwellFD3D"
