module MaxwellFDM

# @reexport makes all exported symbols of the exported packages available in module using MaxwellFDM.
using Reexport
@reexport using LinearAlgebra, SparseArrays, StaggeredGridCalculus, GeometryPrimitives
using StaticArrays

export SComplex, SSComplex3, ParamInd, ObjInd

## Type aliases
# Below, use Int instead of Int64 for compatibility with 32-bit systems (e.g., x86 in appveyor.yml).
const Float = typeof(0.0)  # use Float = Float128 for quadruple precision in the future
const CFloat = Complex{Float}

const Tuple2 = NTuple{2}
const Tuple3 = NTuple{3}
const Tuple4 = NTuple{4}
Tuple22{T} = Tuple2{Tuple2{T}}
Tuple23{T} = Tuple2{Tuple3{T}}
Tuple24{T} = Tuple2{Tuple4{T}}
Tuple32{T} = Tuple3{Tuple2{T}}
Tuple33{T} = Tuple3{Tuple3{T}}

const AbsVec = AbstractVector
const AbsMat = AbstractMatrix
const AbsArr = AbstractArray

const VecBool = Vector{Bool}
const VecInt = Vector{Int}
const VecFloat = Vector{Float}
const VecComplex = Vector{CFloat}

const AbsVecBool = AbsVec{Bool}
const AbsVecInt = AbsVec{Int}
const AbsVecFloat = AbsVec{Float}
const AbsVecComplex = AbsVec{CFloat}

const AbsVecInteger = AbsVec{<:Integer}
const AbsVecReal = AbsVec{<:Real}
const AbsVecNumber = AbsVec{<:Number}

const MatFloat = Matrix{Float}
const MatComplex = Matrix{CFloat}

const AbsMatFloat = AbsMat{Float}
const AbsMatComplex = AbsMat{CFloat}

const AbsMatReal = AbsMat{<:Real}
const AbsMatNumber = AbsMat{<:Number}

const AbsArrComplex = AbsArr{CFloat}
const AbsArrNumber{N} = AbsArr{<:Number,N}

const SBool{K} = SVector{K,Bool}
const SFloat{K} = SVector{K,Float}
const SInt{K} = SVector{K,Int}
const SComplex{K} = SVector{K,CFloat}
const SNumber{K} = SVector{K,<:Number}

const SSComplex{K,L} = SMatrix{K,K,CFloat,L}

const SSComplex1 = SSComplex{1,1}
const SSComplex2 = SSComplex{2,4}
const SSComplex3 = SSComplex{3,9}

const MatParam = Union{Number,AbsVecNumber,AbsMatNumber}
const ParamInd = UInt8  # change this to handle more than 2⁸ = 256 materials
const ObjInd = UInt16  # change this to handle more than 2¹⁶ = 65536 objects


# The order of inclusion matters: if types or functions in file A are used in file B, file A
# must be included first.
include("enumtype.jl")
include("util.jl")
include("phys.jl")
include("material.jl")
include("object.jl")
include("field.jl")
include("source/source.jl")
include("assignment.jl")
include("smoothing.jl")
include("param.jl")
# include("equation.jl")
include("maxwell.jl")

end # module MaxwellFDM
