module MaxwellFD3D

using StaticArrays, GeometryPrimitives
using DataStructures  # SortedSet in gridgen.jl

export SVec3Complex, SMat3Complex, ParamInd, ObjInd

## Type aliases
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
const AbsVecInteger = AbsVec{<:Integer}
const AbsVecReal = AbsVec{<:Real}
const AbsVecNumber = AbsVec{<:Number}
const AbsMatNumber = AbsMat{<:Number}
const SVec3Float = SVector{3,Float}
const SVec2Int = SVector{2,Int}
const SVec3Int = SVector{3,Int}
const SVec3Complex = SVector{3,CFloat}
const SMat3Complex = SMatrix{3,3,CFloat,9}
const MatParam = Union{Number,AbsVecNumber,AbsMatNumber}
const ParamInd = UInt8  # change this to handle more than 2⁸ = 256 materials
const ObjInd = UInt16  # change this to handle more than 2¹⁶ = 65536 objects


# The order of inclusion matters: if types or functions in file A are used in file B, file A
# must be included first.
include("enumtype.jl")
include("base.jl")
include("phys.jl")
include("grid.jl")
include("material.jl")
include("shape.jl")
# include("source.jl")
include("gridgen.jl")
include("smoothing.jl")
# include("field.jl")
include("operator.jl")
include("pml.jl")
# include("equation.jl")

end # module MaxwellFD3D
