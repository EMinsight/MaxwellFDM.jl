using MaxwellFDM
using Test
using Statistics: mean
using LinearAlgebra, SparseArrays

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellFDM" begin

# include("enumtype.jl")
# include("base.jl")
# include("phys.jl")
# include("grid.jl")
# include("object.jl")
# include("gridgen.jl")
# include("material.jl")
# include("smoothing.jl")
# include("differential.jl")
# include("mean.jl")
# include("param.jl")
include("pml.jl")
# include("source.jl")

# end  # @testset "MaxwellFDM"
