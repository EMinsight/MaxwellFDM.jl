export Grid  # types
export isproperebc, lghost  # functions

struct Ghosted{K}
    l::Tuple2{NTuple{K,VecFloat}}  # l[PRIM][k] = primal vertex locations with ghost points in k-direction
    τl::Tuple2{NTuple{K,VecFloat}}  # l[PRIM][k] = primal vertex locations with transformed ghost points in k-direction
    τind::Tuple2{NTuple{K,VecInt}}  # τind[PRIM][k] = indices of Ghosted.l corresponding to transformed points by boundary conditions
    ∆τ::Tuple2{NTuple{K,VecFloat}}  # ∆[PRIM][k] = amount of shift to get points shifted by BLOCH boundary conditions (Doesn't care transformation by reflection)
end

# Consider storing lvxlbounds, which is
# ([(ldual₀,ldual₁), (ldual₁,ldual₂), ..., (ldualₙ₋₁,ldualₙ)],
#  [(lprim₁,lprim₂), ..., (lprimₙ₋₁,lprimₙ), (lprimₙ,lprimₙ₊₁)]
# The purpose of this field is easy retrieval of voxel bounds surrounding a voxel center.
# For example, for a voxel center l[gt][i], the surrounding voxel bounds are lvxlbounds[gt][i],
# regardless of the value of gt.  You don't need to adjust i by ±1.
#
# Another option.  To get the voxel bounds for voxel centers at l[gt], take l[alter(gt)],
# and prepend or append the ghost point appropriately.  Let the resulting array lvxlbounds.
# Then, for a voxel center l[gt][i], the voxel bounds are lvxlbounds[[i,i+1]], regardless of
# the value of gt.
struct Grid{K}
    axis::SVector{K,Axis}  # axes of this grid (one of X̂, Ŷ, Ẑ, the instances of Axis)
    unit::PhysUnit  # set of physical units used on computational grid
    N::SVector{K,Int}  # N[k] = number of grid cells in k-direction
    L::SVector{K,Float}  # L[k] = length of grid (domain) in k-direction
    l::Tuple2{NTuple{K,VecFloat}}  # l[PRIM][k] = primal vertex locations in k-direction
    ∆l::Tuple2{NTuple{K,VecFloat}}  # ∆l[PRIM][k] = (∆l at primal vertices in w) == diff(l[DUAL][k] including ghost point)
    Npml::Tuple2{SVector{K,Int}}  # Npml[NEG][k] = number of cells inside PML at (-) end in k-direction
    Lpml::Tuple2{SVector{K,Float}}  # Lpml[NEG][k] = thickness of PML at (-) end in k-direction
    lpml::Tuple2{SVector{K,Float}}  # lpml[NEG][k] = location of PML interface at (-) end in k-direction
    ebc::SVector{K,EBC}  # ebc[k] = internal boundary condition in k-direction
    σ::Tuple2{NTuple{K,VecBool}}  # false for non-ghost points exactly on symmetric boundary (PPC or PDC)
    boundstype::SVector{K,GridType}  # boundstype[k] = grid type of k-normal boundary planes
    bounds::Tuple2{SVector{K,Float}}  # bounds[NEG][k] = boundary of domain at (-) end in k-direction
    center::SVector{K,Float}  # center of grid without PML
    ghosted::Ghosted{K}  # data related to ghost points
end

# Constructor for 1D grid: arguments don't have to be tuples.
Grid(axis::Axis, unit::PhysUnit, lprim::AbsVecReal, Npml::Tuple2{<:Integer}, ebc::EBC) =
    Grid(SVector(axis), unit, (lprim,), ([Npml[nN]], [Npml[nP]]), SVector(ebc))

# Constructor for 3D grid: axis is always XYZ.
Grid(unit::PhysUnit, lprim::Tuple3{AbsVecReal}, Npml::Tuple2{AbsVecInteger}, ebc::AbsVec{EBC}) =
    Grid(XYZ, unit, lprim, Npml, ebc)

# Constructor taking non-static vectors.
Grid(axis::AbsVec{Axis}, unit::PhysUnit, lprim::NTuple{K,AbsVecReal},
     Npml::Tuple2{AbsVecInteger}, ebc::AbsVec{EBC}) where {K} =
     Grid(SVector{K}(axis), unit, lprim, SVector{K}.(Npml), SVector{K}(ebc))

# Constructor calling the inner constructor.
function Grid(axis::SVector{K,Axis}, unit::PhysUnit, lprim::NTuple{K,AbsVecReal},
              Npml::Tuple2{SVector{K,<:Integer}}, ebc::SVector{K,EBC}) where {K}
    all(issorted.(lprim)) || throw(ArgumentError("all entry vectors of lprim = $(lprim) must be sorted."))

    # For array inputs, create separate copies.
    lprim = deepcopy(lprim)
    Npml = deepcopy(Npml)

    # Determine the grid type of the domain boundary.
    # If ebc[k] = PPC or BLOCH, lprim[k] is from the negative boundary to the positive boundary of the domain.
    # If ebc[k] = PDC, lprim[k] is from outside the negative boundary to outside the positive boundary.
    boundstype = map(e->(e≠PDC ? PRIM : DUAL), ebc)  # SVector{K,GridType}

    ldual = movingavg.(lprim)  # NTuple{K,VecFloat}
    lbound = map((b,p,d)->(b==PRIM ? p : d), boundstype, SVector(lprim), SVector(ldual))  # SVector{K,Vector{<:Real}}

    # Set bounds, L, lpml, Lpml, and center.
    bounds = (getindex.(lbound,1), getindex.(lbound,endof.(lbound)))  # Tuple2{SVector{K,<:Real}}
    L = bounds[nP] - bounds[nN]  # SVector{K,<:Real}

    lpml = (map((v,n)->v[1+n], lbound, Npml[nN]), map((v,n)->v[end-n], lbound, Npml[nP]))  # Tuple2{SVector{K,<:Real}}
    Lpml = (lpml[nN]-bounds[nN], bounds[nP]-lpml[nP])  # Tuple2{SVector{K,<:Real}}

    center = (lpml[nN] + lpml[nP]) / 2  # SVector{K,<:Real}

    # Prepare lprim and ldual to calculate ∆lprim and ∆ldual.
    # Note that ∆lprim is not ∆'s or lprim, but ∆'s defined at lprim locations.
    for k = 1:K
        if boundstype[k] == PRIM  # primal planes define boundaries
            # To calculate ∆'s at lprim, we need one more ldual outside the negative boundary.
            ldual₀ = ebc[k]==BLOCH ? ldual[k][end]-L[k] : 2bounds[nN][k]-ldual[k][1]  # ebc = PPC if ebc ≠ BLOCH
            prepend!(ldual[k], ldual₀)  # length(ldual[k]) = N[k]+1
        else  # dual planes define boundaries
            # lprim[k] is from outside the negative boundary to outside the positive
            # boundary.  ldual[k] that is calculated as movingavg(lprim[k]) is from the
            # negative boundary to the positive boundary.  In this case, ldual[k] has
            # N[k]+1 entries, and ldual[k][1] is the ghost point.
            shift!(lprim[k])  # length(lprim[k]) = N[k]+1
        end
    end

    # Below, (∆lprim, ∆ldual) is not (diff(lprim), diff(ldual)), but swapped.
    ∆lprim = diff.(ldual)  # length(∆lprim[k]) = N[k]
    ∆ldual = diff.(lprim)  # length(∆ldual[k]) = N[k]
    ∆l = (∆lprim, ∆ldual)

    # Set N (number of grid cells along the axis).
    assert(length.(∆l[nPR]) == length.(∆l[nDL]))  # lprim, ldual, ∆lprim, ∆ldual have the same length
    N = SVector(length.(∆l[nPR]))

    # Find the locations of the ghost points transformed into the domain by the boundary
    # conditions.  Note that regardless of the boundary condition, the primal ghost point is
    # on the positive side (so it is lprim[end]) and the dual ghost point is on the negative
    # side (so it is ldual[1]).
    # Make sure these locations are where the objects assigned to the ghost points are found.
    bcp = ebc.==PPC
    bcd = ebc.==PDC
    bcb = ebc.==BLOCH

    # Construct an instance of Ghosted.
    τind = (map(n->collect(1:n+1), N.data), map(n->collect(1:n+1), N.data))  # Tuple23{VecInt}
    τindg_prim = .!bcb.*N + 1 - bcd  # SVec3Int: N+1 for PPC; N for PDC; 1 for BLOCH
    τindg_dual = bcb.*N + 1 + bcp  # SVec3Int: 2 for PPC; 1 for PDC; N+1 for BLOCH
    τindg = (τindg_prim, τindg_dual)  # Tuple2{SVec3Int}

    τl = (deepcopy(lprim), deepcopy(ldual))

    ∆τ = (zeros.((N+1).data), zeros.((N+1).data))  # Tuple23{VecFloat}

    for k = 1:K
        τind[nPR][k][end] = τindg[nPR][k]
        τind[nDL][k][1] = τindg[nDL][k]

        τl[nPR][k][end] = getindex(lprim[k], τind[nPR][k][end])
        τl[nDL][k][1] = getindex(ldual[k], τind[nDL][k][1])

        ∆τ[nPR][k][end] = -L[k] * bcb[k]
        ∆τ[nDL][k][1] = L[k] * bcb[k]
    end

    ghosted = Ghosted{K}((deepcopy(lprim), deepcopy(ldual)), τl, τind, ∆τ)

    # Construct σ, which is false on symmetry boundary (non-BLOCH boundaries).  Note that
    # σ's are vectors of length(N) and therefore do not include ghost points.
    σprim = ones.(Bool, N.data)
    σdual = ones.(Bool, N.data)
    for k = 1:K
        # Primal grid points are on the symmetry boundary if they are the first (non-ghost)
        # points and PPC boundary is used.
        σprim[k][1] = !bcp[k]

        # Dual grid points are on the symmetry boundary if they are the last (non-ghost)
        # points and PDC boundary is used.
        σdual[k][end] = !bcd[k]
    end
    σ = (σprim, σdual)

    # Eliminate ghost points.  Regardless of the boundary type, the primal ghost point
    # is on the positive side (so lprim[end] is dropped) and the dual ghost point is on
    # the negative side (so ldual[1] is dropped).
    pop!.(lprim)
    shift!.(ldual)
    assert(length.(lprim) == length.(ldual) == N.data)  # lprim, ldual, ∆lprim, ∆ldual have the same length
    l = (lprim, ldual)

    return Grid{K}(axis, unit, N, L, l, ∆l, Npml, Lpml, lpml, ebc, σ, boundstype, bounds, center, ghosted)
end

# To do
# When needed, implement `section` (or `cross_section`) that returns Grid{K₁} from Grid{K₀}
# for K₁ < K₀.  Provide the directions (Axis) across which the cross section in taken.  Use
# the inner constructor of Grid to construct the cross section.

# Additional functions
Base.ndims(::Grid{K}) where {K} = K
Base.in(l::SVector{K}, g::Grid{K}) where {K} = all(g.bounds[nN] .≤ l .≤ g.bounds[nP])

# Check if two ends are both BLOCH or both non-BLOCH.  If both are BLOCH, check if e⁻ⁱᵏᴸ is
# a phase factor (i.e., amplitude = 1).  If both are non-BLOCH, check if e⁻ⁱᵏᴸ = 1.
#
# e⁻ⁱᵏᴸ is the accumulated phase from the negative-end boundary to the positive-end one.  L
# is the distance between the boundaries (= domain size).  It is exp(-ikL), not exp(+ikL),
# because we assume the exp(+iωt) time dependence and thus the position dependence is
# exp(-ik⋅r).
isproperebc(ebc::EBC, e⁻ⁱᵏᴸ::Number=1) = (ebc==BLOCH && abs(e⁻ⁱᵏᴸ)==1) || (ebc≠BLOCH && e⁻ⁱᵏᴸ==1)
isproperebc(ebc::SVector{K,EBC}, e⁻ⁱᵏᴸ::SVector{K,<:Number}=@SVector(ones(K))) where {K} =
    all(ntuple(k->isproperebc(ebc[k], e⁻ⁱᵏᴸ[k]), Val{K}))
