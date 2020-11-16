export create_field3d, field3d2vec

# About the order of indices of f3d:
#
# Like param3d in assignment.jl, we index f3d as f3d[i,j,k,w], where (i,j,k) are positional
# indices and w is the Cartesian component index.
#
# The reason for this choice is the same as param3d.  For example, in assigning source
# values to j3d (an array of current density), we usually fix a component first and then
# assign the same value to a range of (i,j,k).  This can be more efficiently done with
# j3d[i,j,k,w], because a contiguous range of (i,j,k) actually corresponds to a contiguous
# memory block.
#
# I think I will need to index the E- and H-fields the same way for matrix-free operations.
# When we perform curl operation on these field arrays, we implement the ∂/∂w operation as
# fixing the component to differentiate first and then perform the differentiation.
# Therefore, again the E[i,j,k,w] indexing scheme results in an operation on a more
# contiguous block in memory space.

create_field3d(N::SInt{3}) = zeros(CFloat, N.data..., 3)  # 3 = numel(Axis)

# Below, permutedims(f3d, ...) create a new array, whereas reshape(f3d, :) doesn't.
# Therefore, if implemented naively, this function creates a new array for order_cmpfirst =
# true whereas it doesn't for order_cmpfirst = false.
field3d2vec(f3d::AbsArrNumber{4}; order_cmpfirst::Bool=true) =
    order_cmpfirst ? reshape(permutedims(f3d, (4,1,2,3)), :) : reshape(f3d,:)
