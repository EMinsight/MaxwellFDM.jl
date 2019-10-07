# import Base:isapprox, dot

# Below, t_ind returns
# - a tuple if the first argument is a tuple of tuple, and
# - an SVector if the first argument is a tuple of vector.

# Below, earlier methods delegate actions to later methods.  Maybe they need to be
# implemented separately for speed?

# t_ind(t::Tuple3{T}, ind::Tuple3{Int}) where {T} = (t[ind[1]], t[ind[2]], t[ind[3]])
# t_ind(t::Tuple32{T}, i::Int, j::Int, k::Int) where {T} = (t[1][i], t[2][j], t[3][k])  # i, j, k = 1 or 2

# From a tuple of two tuples 1 and 2, each with K entries, construct a tuple (with K entries)
# whose kth entry is the kth entry of either tuple 1 or tuple 2.
# E.g., t_ind(((0.1,0.2,0.3), (1.0,2.0,3.0)), 1, 2, 1) = (0.1, 2.0, 0.3)
@inline t_ind(t::Tuple23, i₁₂::T, j₁₂::T, k₁₂::T) where {T<:Union{GridType,Sign,Integer}} = t_ind(t, (i₁₂,j₁₂,k₁₂))  # i₁₂, j₁₂, k₁₂ = 1 or 2
@inline t_ind(t::Tuple2{NTuple{K}}, ind₁₂::CartesianIndex{K}) where {K} = t_ind(t, ind₁₂.I)  # ind₁₂.I: NTuple{K,Int}
@inline t_ind(t::Tuple2{NTuple{K}}, ind₁₂::SVector{K,T}) where {K,T<:Union{GridType,Sign,Integer}} = t_ind(t, ind₁₂.data)
@inline t_ind(t::Tuple2{NTuple{K}}, ind₁₂::NTuple{K,T}) where {K,T<:Union{GridType,Sign,Integer}} =  # ind₁₂[k] = 1 or 2
    map((t₁,t₂,i) -> Int(i)==1 ? t₁ : t₂, t[1], t[2], ind₁₂)  # NTuple{K}

# From a tuple of two SVectors 1 and 2, each with K entries, construct an SVector (with K
# entries) whose kth entry is the kth entry of either SVector 1 or SVector 2.
# E.g., t_ind((SVector(0.1,0.2,0.3), SVector(1.0,2.0,3.0)), 1, 2, 1) = SVector(0.1, 2.0, 0.3)
@inline t_ind(t::Tuple2{SVec3}, i₁₂::T, j₁₂::T, k₁₂::T) where {T<:Union{GridType,Sign,Integer}} = t_ind(t, (i₁₂,j₁₂,k₁₂))
@inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::CartesianIndex{K}) where {K} = t_ind(t, ind₁₂.I)
@inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::NTuple{K,T}) where {K,T<:Union{GridType,Sign,Integer}} = t_ind(t, SVector(ind₁₂))
@inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::SVector{K,T}) where {K,T<:Union{GridType,Sign,Integer}} =  # ind₁₂[k] = 1 or 2
    map((t₁,t₂,i) -> Int(i)==1 ? t₁ : t₂, t[1], t[2], ind₁₂)  # SVector{K}

# From a tuple of K vectors, construct a vector with K entries whose kth entry is
# taken from the kth vector.  Which entry to take from the kth vector is specified
# by indices.
# E.g., t_ind(([0.1,0.2,0.3,0.4], [1.0,2.0,3.0,4.0], [10.0,20.0,30.0,40.0]), 3, 1, 4) = SVector(0.3, 1.0, 40.0)
@inline t_ind(t::Tuple3{AbsVec}, i::Int, j::Int, k::Int) = t_ind(t, (i,j,k))
@inline t_ind(t::Tuple2{AbsVec}, i::Int, j::Int) = t_ind(t, (i,j))
@inline t_ind(t::NTuple{K,AbsVec}, ind::CartesianIndex{K}) where {K} = t_ind(t, ind.I)
@inline t_ind(t::NTuple{K,AbsVec}, ind::NTuple{K,Int}) where {K} = t_ind(t, SVector(ind))
@inline t_ind(t::NTuple{K,AbsVec}, ind::SVector{K,Int}) where {K} = map((tₖ,iₖ) -> tₖ[iₖ], SVector(t), ind)

# getindex(t::Tuple3{T}, ind::Tuple3{Int}) where {T} = (t[ind[1]], t[ind[2]], t[ind[3]])
# getindex(t::Tuple32{T}, i::Int, j::Int, k::Int) where {T} = (t[1][i], t[2][j], t[3][k])  # i, j, k = 1 or 2
# getindex(t::Tuple23{T}, i::Int, j::Int, k::Int) where {T} = (t[i][1], t[j][2], t[k][3])  # i, j, k = 1 or 2
# getindex(t::Tuple3{AbsVec{T}}, i::Int, j::Int, k::Int) where {T} = (t[1][i], t[2][j], t[3][k])
# getindex(t::Tuple3{AbsVec{T}}, i::Tuple2{Int}, j::Tuple2{Int}, k::Tuple2{Int}) where {T} =
#     ((t[1][i[1]], t[1][i[2]]), (t[2][j[1]], t[2][j[2]]), (t[3][k[1]], t[3][k[2]]))

# rand3() = (rand(), rand(), rand())
# randn3() = (randn(), randn(), randn())
#
# dot(t1::Tuple3{Real}, t2::Tuple3{Real}) = t1[1]*t2[1] + t1[2]*t2[2] + t1[3]*t2[3]
#
# function norm2(x)
#     if isa(x, Number)
#         return abs(x)
#     else
#         sx = start(x)
#
#         result = 0.0
#         (xi, sx) = next(x, sx)
#         result += norm2(xi)^2
#         while !done(x, sx)
#             (xi, sx) = next(x, sx)
#             result += norm2(xi)^2
#         end
#
#         return sqrt(result)
#     end
# end
#
# function norm2diff(x, y)
#     if isa(x, Number)
#         if !isa(y, Number)
#             return Inf
#         else
#             return abs(x-y)
#         end
#     elseif (isa(x,AbsArr) && isa(y,Tuple)) || (isa(x,Tuple) && isa(y,AbsArr))
#         return Inf
#     else
#         sx, sy = start(x), start(y)
#
#         result = 0.0
#         (xi, sx), (yi, sy) = next(x, sx), next(y, sy)
#         result += norm2diff(xi, yi)^2
#         while !done(x, sx) && !done(y, sy)
#             (xi, sx), (yi, sy) = next(x, sx), next(y, sy)
#             result += norm2diff(xi, yi)^2
#         end
#
#         if !done(x, sx) || !done(y, sy)
#             return Inf
#         else
#             return sqrt(result)
#         end
#     end
# end
