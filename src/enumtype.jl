export FieldType, BC
export EH, nE, nH, nEH

# Field types
const nE, nH = 1, 2  # E-, H-fields
const nEH = SVector(nE, nH)
@enum FieldType EE=nE HH
const EH = SVector(EE, HH)
for ins in instances(FieldType); @eval export $(Symbol(ins)); end  # export all instances
Base.string(ins::FieldType) = ins==EE ? "E" : "H"
StaggeredGridCalculus.alter(ins::FieldType) = ins==EE ? HH : EE

# Given boundary field types, determine whether the grid planes specified by the given field
# type ft are primal or dual grid planes in the Cartesian directions.
ft2gt(ft::FieldType, boundft::FieldType) = PD[2 - (boundft==ft)]

# Return grid types of the w-component of the field F.
function gt_Fw(nw::Int,  # component index of field
               ft::FieldType,  # type of field
               boundft::SVector{K,FieldType}) where {K}  # type of field parallel to domain boundaries
    gt₀ = ft2gt.(ft, boundft)  # grid type of corners of voxel whose edges are ft lines
    gt = broadcast((a,w,g)->(a==w ? alter(g) : g), SVector(ntuple(identity, Val(K))), nw, gt₀)  # grid type of Fw; no change from gt₀ for w ∉ axes

    return gt
end

# Boundary conditions
@enum BC PERIODIC=1 CONDUCTING  # periodic boundary condition, perfect conductor boundary condition
for ins in instances(BC); @eval export $(Symbol(ins)); end  # export all instances
Base.string(ins::BC) = ins==PERIODIC ? "periodic" : "perfect conductor"
