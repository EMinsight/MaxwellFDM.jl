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
@inline t_ind(t::Tuple2{SVector{3}}, i₁₂::T, j₁₂::T, k₁₂::T) where {T<:Union{GridType,Sign}} = t_ind(t, SVector(i₁₂,j₁₂,k₁₂))
# @inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::CartesianIndex{K}) where {K} = t_ind(t, ind₁₂.I)
@inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::NTuple{K,T}) where {K,T<:Union{GridType,Sign}} = t_ind(t, SVector(ind₁₂))
@inline t_ind(t::Tuple2{SVector{K}}, ind₁₂::SVector{K,T}) where {K,T<:Union{GridType,Sign}} =  # ind₁₂[k] = 1 or 2
    map((t₁,t₂,i) -> Int(i)==1 ? t₁ : t₂, t[1], t[2], ind₁₂)  # SVector{K}

# From a tuple of K vectors, construct a vector with K entries whose kth entry is
# taken from the kth vector.  Which entry to take from the kth vector is specified
# by indices.
# E.g., t_ind(([0.1,0.2,0.3,0.4], [1.0,2.0,3.0,4.0], [10.0,20.0,30.0,40.0]), 3, 1, 4) = SVector(0.3, 1.0, 40.0)
@inline t_ind(t::Tuple3{AbsVec}, i::Int, j::Int, k::Int) = t_ind(t, SVector(i,j,k))
@inline t_ind(t::Tuple2{AbsVec}, i::Int, j::Int) = t_ind(t, SVector(i,j))
@inline t_ind(t::NTuple{K,AbsVec}, ind::CartesianIndex{K}) where {K} = t_ind(t, ind.I)
@inline t_ind(t::NTuple{K,AbsVec}, ind::NTuple{K,Int}) where {K} = t_ind(t, SVector(ind))
@inline t_ind(t::NTuple{K,AbsVec}, ind::SVector{K,Int}) where {K} = map((tₖ,iₖ) -> tₖ[iₖ], SVector(t), ind)  # SVector{K}

# See http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array.
# Generate the swap macros at http://pages.ripco.net/~jgamble/nw.html.
# - Number of inputs: length of the input vector
# - Algorithm choices: "Best"
# - Select "Create a set of SWAP macros".
@inline function swap_ind!(ind::AbsArrInteger, v::AbsArr, i::Integer, j::Integer)
    @inbounds if v[ind[j]] < v[ind[i]]
        @inbounds ind[i], ind[j] = ind[j], ind[i]
    end
end

@inline function sort_ind!(ind::MArray{Tuple{2,2,2},<:Integer,3,8}, v::MArray{Tuple{2,2,2},<:Real,3,8})
    # Initialize ind.
    @simd for n = 1:8
        @inbounds ind[n] = n
    end

    # Sort ind.
    swap_ind!(ind,v,1,2); swap_ind!(ind,v,3,4); swap_ind!(ind,v,1,3); swap_ind!(ind,v,2,4);
    swap_ind!(ind,v,2,3); swap_ind!(ind,v,5,6); swap_ind!(ind,v,7,8); swap_ind!(ind,v,5,7);
    swap_ind!(ind,v,6,8); swap_ind!(ind,v,6,7); swap_ind!(ind,v,1,5); swap_ind!(ind,v,2,6);
    swap_ind!(ind,v,2,5); swap_ind!(ind,v,3,7); swap_ind!(ind,v,4,8); swap_ind!(ind,v,4,7);
    swap_ind!(ind,v,3,5); swap_ind!(ind,v,4,6); swap_ind!(ind,v,4,5)

    return ind
end

@inline function sort_ind!(ind::MArray{Tuple{2,2},<:Integer,2,4}, v::MArray{Tuple{2,2},<:Real,2,4})
    # Initialize ind.
    @simd for n = 1:4
        @inbounds ind[n] = n
    end

    # Sort ind.
    swap_ind!(ind,v,1,2); swap_ind!(ind,v,3,4); swap_ind!(ind,v,1,3); swap_ind!(ind,v,2,4); swap_ind!(ind,v,2,3);

    return ind
end

@inline function sort_ind!(ind::MArray{Tuple{2},<:Integer,1,2}, v::MArray{Tuple{2},<:Real,1,2})
    # Initialize ind.
    @simd for n = 1:2
        @inbounds ind[n] = n
    end

    # Sort ind.
    swap_ind!(ind,v,1,2);

    return ind
end

isuniform(mv::MArray{Tuple{2,2,2},<:Any,3,8}) = @inbounds (mv[1]==mv[2]==mv[3]==mv[4]==mv[5]==mv[6]==mv[7]==mv[8])
isuniform(mv::MArray{Tuple{2,2},<:Any,2,4}) = @inbounds (mv[1]==mv[2]==mv[3]==mv[4])
isuniform(mv::MArray{Tuple{2},<:Any,1,2}) = @inbounds (mv[1]==mv[2])

# Determine if the field subspace is orthogonal to the shape subspace; used in smoothing.jl
# and source/*.jl.
#
# The shape and field subspaces are the subspaces of the 3D space where the shapes and
# fields lie.  In standard 3D problems, both subspaces are 3D.  However there are other
# cases as well: for example, in 2D TM problems, the shapes are on the xy-plane (2D space),
# but the E-fields are along the z-direction (1D space).
#
# For the shape and field subspace dimensions K and Kf, orthogonal subspaces can occur only
# when K + Kf ≤ 3, because if two subspaces are orthogonal to each other and K + Kf > 3, we
# can choose K + Kf linearly independent vectors in the direct sum of the two subspaces,
# which is impossible as the direct sum should still be a subspace of the 3D space.
#
# Considering the above, there are only a few combinations of K and Kf that allow orthogonal
# subspaces: (K,Kf) = (2,1), (1,1), (1,2).
# - (K,Kf) = (2,1).  This happens in 2D TE or TM problems.  The shape subspace is the 2D xy-
# plane, but the magnetic (electric) field subspace in TE (TM) problems is the 1D z-axis.
# - (K,Kf) = (1,1).  This happens in the 1D TEM problems with isotropic materials.  The
# shape subspace is the 1D z-axis, but the E- and H-field spaces are the 1D x- and y-axes.
# - (K,Kf) = (1,2).  This happens in the 1D TEM problem with anisotropic materials.
#
# Note that we can always solve problems as if they are 3D problems.  So, the consideration
# of the cases with K, Kf ≠ 3 occurs only when we can build special equations with reduced
# number of degrees of freedom, like in the TE, TM, TEM equations.  In such cases, we find
# that the two subspaces are always orthogonal if K ≠ Kf.  In fact, smooth_param!() is
# written such that the Kottke's subpixel smoothing algorithm that decomposes the field into
# the components tangential and normal to the shape surface is applied only when the shape
# and field subspaces coincide.  (This makes sense, because the inner product between the
# field and the direction normal that needs to be performed to decompose the field into such
# tangential and normal components is defined only when the field and the direction normal
# are in the same vector space.)  Therefore, if we want to apply Kottke's subpixel smoothing
# algorithm that depends on the surface normal direction, we have to construct the problem
# such that K = Kf (but the converse is not true: K = Kf does not imply the use of Kottke's
# subpixel smoothing algorithm that depends on the surface normal direction; in other words,
# the case with K = Kf can still use the subpixel smoothing algorithm that assumes the
# orthogonality between the shape and field subspace, such that the field is always
# tangential to the shape surface.)
#
# The contraposition of the above statement is that if K ≠ Kf, then Kottke's subpixel
# smoothing algorithm that does NOT depend on the direction normal is applied.  This
# subpixel smoothing algorithm assumes that the field subspace is orthogonal to the shape
# subspace, such that the field has only the tangential component to the surface.
# Therefore, if we pass K ≠ Kf to smooth_param!(), we should make sure that the field
# subspace is orthogonal to the shape subspace.  Note that this does not exclude the
# possibility of K ≠ Kf while the two subspaces are nonorthogonal; it just means that we
# must formulate the problem differently in such cases by, e.g., decomposing the equations,
# in order to use smooth_param!().
#
# As noted earlier, K = Kf can stil include cases where the shape and field subspaces are
# orthogonal.  Because K + Kf = 2K should be ≤ 3 in+ such cases, we conclude that K = Kf = 1
# is the only case where the two subspaces could be orthogonal while K = Kf.  In fact, when
# K = Kf, we will assume that the two subspaces are orthogonal, because the problems with 1D
# slabs and the field along the slab thickness direction are not interesting from the EM
# wave propagation perspectives.
isfield_ortho_shape(Kf, K) = Kf≠K || Kf==1
