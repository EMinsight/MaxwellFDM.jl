using DataStructures: SortedSet

export gen_lprim
export R_TARGET_DEFAULT, R_MAX_DEFAULT

const R_TARGET_DEFAULT = 1.9  # default target ratio of geometric sequence of grid cell size; use 1.2 to be more conservative
const R_MAX_DEFAULT = 2.0  # maximum ratio of geometric sequence of grid cell size; use 1.3 to be more conservative


function movingmin(l::AbsVec{T}) where {T<:Number}
    # Return min.([l; Inf], [Inf; l])

    n = length(l)
    if n == 0
        lmov = [Inf]
    else  # n ≥ 1
        lmov = Vector{T}(undef, n+1)
        lmov[1] = l[1]
        lmov[n+1] = l[n]
        for i = 2:n
            lmov[i] = min(l[i-1], l[i])
        end
    end

    return lmov
end

function consolidate_similar!(l::AbsVecNumber, L::Number)
    # Consolidate approximately equal points.  Assume l is sorted.

    # Original code:
    # ind_unique = (≉).(@view(lprim0[1:end-1]), @view(lprim0[2:end]))  # length(ind_unique) = length(lprim0) - 1
    # push!(ind_unique, true)
    # lprim0 = lprim0[ind_unique]  # [lprim0[1:end-1][original ind_unique]; lprim0[end]]

    if !isempty(l)
        ind_sim = findall(isapprox_wrt.(@view(l[1:end-1]), @view(l[2:end]), L))

        # Below, l[end] must be always included as a unique entry.  There are two cases: where
        # l[end-1] is included as a unique entry or not.  If l[end-1] is included, it must be
        # very different from l[end], so we must include l[end].  If l[end-1] is not included,
        # it must be approximately l[end], but because l[end-1] is not included, we must include
        # l[end].  In other words, l[end] represents the last set of approximately equal entries
        # in any cases.
        deleteat!(l, ind_sim)
        @assert !isapprox_wrt(l[end], l[end-1], L)  # last entry is always unique
    end
end

function findsim(l1::AbsVecNumber, l2::AbsVecNumber, L::Number)
    ind = Int[]
    for i = 1:length(l1)
        v1 = l1[i]
        for v2 = l2
            isapprox_wrt(v1, v2, L) && push!(ind, i)
        end
    end

    return ind
end

# From a given array of ClosedInterval's, generate a subset of primal grid.
function gen_sublprim1d(domain::OpenInterval,  # specifies domain boundaries; used to ignore points outside domain
                        domaintype::GridType,
                        intvprim::AbsVec{ClosedInterval},  # object boundaries to fit to primal grid planes; does not change
                        intvdual::AbsVec{ClosedInterval},  # object boundaries to fit to dual grid planes; does not change
                        lprim₀::AbsVecFloat,  # default points to take as primal grid planes; copied internally
                        ldual₀::AbsVecFloat)  # default points to take as dual grid planes; copied internally

    # (length(Lpml)==2 && all(Lpml.≥0)) || throw(ArgumentError("Lpml = $Lpml must be length-2 and nonnegative."))
    # sum(Lpml) ≤ L_(domain) || throw(ArgumentError("sum(Lpml) = $(sum(Lpml)) must be ≤ L_(domain) = $(L_(domain))."))

    # The following is not necessary, because we will simply ignore points outside the domain.
    # all(contains.(domain, lprim0)) || throw(ArgumentError("Entries of lprim0 = $lprim0 must be all inside domain = $domain."))
    # all(contains.(domain, ldual0)) || throw(ArgumentError("Entries of ldual0 = $ldual0 must be all inside domain = $domain."))

    ∆lmax = max∆l(domain)
    L = length(domain)

    # Collect boundary points.
    # Lpml = float(Lpml)
    # Lpml_signed = [Lpml[1], -Lpml[2]]
    # lprim0 = [lprim0; domain.bound...; [domain.bound...]+Lpml_signed]  # possible duplicates in case of Lpml == 0 are removed later

    # The line below is commented out, because the domain boundaries can be dual grid locations for the PDC boundary condition.
    # lprim0 = [lprim0; domain.bound...]  # possible duplicates in case of Lpml == 0 are removed later

    # Below, exclude the domain boundaries from testing containment, because they will be added later.
    lprim0set = SortedSet(lprim₀)
    for intv = intvprim
        b = bounds(intv)
        b[1] ∈ domain && push!(lprim0set, b[1])
        b[2] ∈ domain && push!(lprim0set, b[2])
    end
    if domaintype == PRIM
        b = bounds(domain)
        push!(lprim0set, b[1], b[2])
    end
    lprim0 = collect(lprim0set)
    consolidate_similar!(lprim0, L)

    # Repeat for dual.
    ldual0set = SortedSet(ldual₀)
    for intv = intvdual
        b = bounds(intv)
        b[1] ∈ domain && push!(ldual0set, b[1])
        b[2] ∈ domain && push!(ldual0set, b[2])
    end
    if domaintype == DUAL
        b = bounds(domain)
        push!(ldual0set, b[1], b[2])
    end
    ldual0 = collect(ldual0set)
    consolidate_similar!(ldual0, L)

    indsim = findsim(ldual0, lprim0, L)
    !isempty(indsim) && @warn "$(ldual0[indsim]) is shared by primal and dual grid, and will be assigned to primal grid."
    deleteat!(ldual0, indsim)

    # Below, insert primal nodes around dual grid points.  This is necessary, because
    # eventually we want to generate only the primal grid; the dual grid is simply
    # defined by the midpoints between the primal nodes.  However, there are some
    # dual nodes that we know must exist in the dual grid.  Therefore, we must
    # create the primal grid such that the dual grid it induces contains those nodes.

    # Make ∆ldual0[j] the smaller of the two ∆ldual0 entries adjacent to ldual0[j].
    ∆ldual0 = movingmin(diff(ldual0))

    intervals = [intvprim; intvdual]
    nprim0 = length(lprim0)
    lprimset_about_ldual0 = SortedSet(Float[])
    for j in eachindex(ldual0)
        l::Float = ldual0[j]
        ind₋ = findlast(x->x<l, lprim0)  # to find interval ∆lprim0[ind] containing ldual0[j]
        ind₋==nothing && (ind₋ = 0)

        # Below, ind₊ is findfirst(x->x>l, lprim0), but because l cannot be one of the
        # lprim0 entries, it can be more efficiently calculated.
        ind₊ = ind₋==nprim0 ? 0 : ind₋+1

        # Find the minimax (minimum of maximum) ∆lprim's (constrained by lprim0 and ldual0) to use around ldual0[j].
        if ind₋ == ind₊ == 0
            @assert nprim0 == 0
            ∆lmm = min(∆lmax, ∆ldual0[j])
        elseif ind₋==0 && ind₊≠0  # ind₊ = 1
            ∆lmm = min(∆lmax, ∆ldual0[j], 2(lprim0[ind₊]-l))
        elseif ind₋≠0 && ind₊==0  # ind₋ = nprim0
            ∆lmm = min(∆lmax, ∆ldual0[j], 2(l-lprim0[ind₋]))
        else  # ind₋≠0 && ind₊≠0
            ∆lmm = min(∆lmax, ∆ldual0[j], 2(l-lprim0[ind₋]), 2(lprim0[ind₊]-l))
        end

        # ClosedInterval containing ldual0[j] reduce the minimum ∆lmm to use around ldual0[j]
        # even further by their ∆lmax constraints.
        for intv = intervals
            if l ∈ intv
                ∆l = max∆l(intv)
                ∆l < ∆lmm && (∆lmm = ∆l)
            end
        end

        # Using the found ∆lmm, put the primal nodes around ldual0[j], such that ldual0[j]
        # is generated as the midpoint between those primal nodes.

        # If domaintype = PRIM, the dual point l is always strictly inside the domain
        # (excluding the boundary), and ∆lmm is such that l ± ∆lmm/2 is within the domain
        # (including the boundary).  Therefore, l ± ∆lmm/2 must be always added as new primal
        # points.
        #
        # If domaintype = DUAL, the smallest and largest dual points l are on the boundary
        # of the domain, and l ± ∆lmm/2 include ghost points, and not beyond them.  Therefore,
        # Again l ± ∆lmm/2 must be always added as new primal points.
        push!(lprimset_about_ldual0, l - ∆lmm/2, l + ∆lmm/2)
    end

    lprim0set = SortedSet(lprim0)
    union!(lprim0set, lprimset_about_ldual0)
    lprim0 = collect(lprim0set)
    consolidate_similar!(lprim0, L)

    # For each space between primal nodes, find the smallest ∆l suggested by intervals.
    nprim0 = length(lprim0)
    btwn_lprim0 = VecFloat(undef, nprim0-1)  # midpoints between primal nodes
    ∆lprim0 = VecFloat(undef, nprim0-1)
    for i = 1:nprim0-1
        btwn_lprim0[i] = (lprim0[i] + lprim0[i+1]) / 2
        ∆lprim0[i] = lprim0[i+1] - lprim0[i]
    end
    # btwn_lprim0 = (lprim0[1:end-1] + lprim0[2:end]) / 2  # points between primal nodes
    # ∆lprim0 = diff(lprim0)

    ∆lmm_array = VecFloat(undef, nprim0-1)  # lmm for each space between primal nodes
    for j = 1:nprim0-1
        l = btwn_lprim0[j]
        ∆lmm = min(∆lmax, ∆lprim0[j])  # minimax: minimum of maximum ∆lprim's allowed

        # l is a midpoint between two primal points, and it can be a dual point by
        # construction of primal points.  Below, the assumption is that such a dual point is
        # contained in the dual intervals (as a boundary point of an interval) in intvdual,
        # such that ∆lmm is correctly shrunken by the ∆lmax constraints of the dual intervals.
        for intv = intervals
            if l ∈ intv
                ∆l = max∆l(intv)
                ∆l < ∆lmm && (∆lmm = ∆l)
            end
        end

        ∆lmm_array[j] = ∆lmm
    end

    # ∆lmm_array is ∆l in the space between primal planes.  From ∆lmm_array, determine ∆l at
    # primal nodes (i.e., a single ∆l to use on the left and right of a primal node).
    ∆l_at_lprim0 = movingmin(∆lmm_array)
    # ∆l_at_lprim0 = min.([Inf; ∆lmm_array], [∆lmm_array; Inf])

    # Create subgrids.  (This is the part that needs to be improved.)
    l = lprim0[1]
    # @assert l == domain.bound[1]  # legacy code before using dual grid for PDC boundary

    ∆l = ∆l_at_lprim0[1]
    prev = [l, l+∆l]
    sublprim = [prev]
    for j = 2:nprim0
        l = lprim0[j]
        ∆l = ∆l_at_lprim0[j]
        if j == nprim0
            curr = [l-∆l, l]
        else  # l is not the first or last lprim0 entries
            curr = [l-∆l, l, l+∆l]  #  same cell size on both sides of primal node; important for averaging ε simply as ε = (ε₁+ε₂)/2
        end

        if isapprox_wrt(@view(curr[1:2]), @view(prev[end-1:end]), L)  # not element-wise comparison
            append!(prev, @view(curr[3:end]))  # exclude curr[1:2]
            # curr = [prev; curr[3:end]]  # exclude curr[1:2]
            # sublprim[end] = curr
        elseif isapprox_wrt(curr[1], prev[end], L)
            append!(prev, @view(curr[2:end]))  # exclude curr[1]
            # curr = [prev; curr[2:end]]  # exclude curr[1]
            # sublprim[end] = curr
        else
            # The followings are equivalent to sublprim = [sublprim, [∆lmm_array[[j-1]], curr]]
            push!(sublprim, [∆lmm_array[j-1]])  # ∆lmm_array[j-1] is scalar; [∆lmm_array[j-1]] is vector with 1 entry
            push!(sublprim, curr)
            prev = curr
        end
    end

    return sublprim
end

function findfirst_stiff_∆∆l(∆l::AbsVecReal, rt::Real)
    n = length(∆l)
    n ≥ 2 || throw(ArgumentError("∆l = $∆l must be length-2 or longer."))

    rt ≥ 1.0 || (rt = 1/rt)  # make rt greater than 1

    for ind = 1:n-1
        (∆l[ind]/∆l[ind+1] > rt || ∆l[ind]/∆l[ind+1] < 1/rt) && return ind
    end

    return 0
end

function issmooth(∆l::AbsVecReal, rt::Real)
    return findfirst_stiff_∆∆l(∆l, rt) == 0  # no stiff_∆∆l
end

function fill_constant(∆lout::Real, gap::AbsVecReal, ∆lt::Real, rt::Real=R_TARGET_DEFAULT, rmax::Real=R_MAX_DEFAULT)
    # Fill the gap with constant ∆l close to ∆lt.

    # ∆lout: ∆l outside gap
    # ∆lt: target constant ∆l inside gap

    ∆lout ≤ ∆lt || ∆lout ≈ ∆lt || throw(ArgumentError("∆lout = $∆lout must be less than or equal to ∆lt = $∆lt."))
    issmooth(SVector(∆lout, ∆lt), rt) || throw(ArgumentError("∆lout = $∆lout and ∆lt = $∆lt must not be too different."))
    rt ≤ rmax || throw(ArgumentError("rt = $rt must be less than or equal to rmax = $rmax"))

    length(gap) == 2 || throw(ArgumentError("gap = $gap must be length-2."))
    (L = diff(gap)[1]) > 0 || throw(ArgumentError("gap[2] = $(gap[2]) must be strictly greater than gap[1] = $(gap[1])."))

    # Below, by taking [floor(), ceil()] instead of [floor(), floor()+1] or [ceil()-1, ceil()],
    # we ensure only one ∆l is considered when L is an exact multiple of ∆lt.
    ns = SVector(floor(Int, L/∆lt), ceil(Int, L/∆lt))  # candidate numbers of grid cells
    ∆l = L ./ ns  # ∆l[1] ≈ ∆l[2] ≈ ∆lt, but ∆l[1] is coarser grid spacing than ∆l[2]

    # ∆lt must be a smooth change from ∆lout (with a rate of change ≤ rt).  We construct
    # ∆l[1] and ∆l[2] close to ∆lt, but neither of them is likely to be exactly ∆lt and
    # therefore could be a slightly stiff change from ∆lout (with a rate of change slightly
    # greater than rt).  Still, if the rate of change is ≤ rmax, we consider ∆l[1] or ∆l[2]
    # acceptable.  If both ∆l[1] and ∆l[2] has a rate of change > rmax from ∆lout, generate
    # a warning.
    issmooth(SVector(∆lout, ∆l[1]), rmax) || issmooth(SVector(∆lout, ∆l[2]), rmax) ||
        # throw(ArgumentError("Cannot find ∆l whose multiple fits in gap = $gap while being close to ∆lt = $∆lt and varying smoothly from ∆lout = $∆lout."))
        @warn "Cannot find ∆l whose multiple fits in gap = $gap while being close to ∆lt = $∆lt and varying smoothly from ∆lout = $∆lout."

    # Unless the coarser grid fails to generate smoothly varying grid, give it a preference.
    i = issmooth(SVector(∆lout, ∆l[1]), rmax) ? 1 : 2
    n = ns[i]
    # ∆l[i] ≤ ∆lt || @info "gap = $gap will be filled with ∆l = $(∆l[i]), which is slightly greater than ∆lt = $∆lt."

    filler = range(gap[1], stop=gap[2], length=n+1)  # in exact arithmetic, equivalent to gap[1]:∆l:gap[end], which may not include gap[end] due to round-off error

    return filler
end

function fill_geometric_sym(∆lsym::Real, gap::AbsVecReal, ∆lt::Real, rt::Real=R_TARGET_DEFAULT, rmax::Real=R_MAX_DEFAULT)
    # Generate a filler subgrid whose grid node separations increases geometrically
    # from both ends.  The target ∆l must be greater than or equal to ∆ln, which is before
    # the gap, and ∆lp, which is after the gap, and ∆ln == ∆lp.

    # gap: length-2 row vector whose elements are the locations of the beginning and
    # ending grid vertices of the gap
    # ∆lsym: ∆l of the subgrids before and after the gap, i.e., ∆ln == ∆lp
    # ∆lt: target ∆l in the gap
    # rt: target ratio of geometric sequence
    # rmax: maximum ratio of geometric sequence

    ∆lt ≥ ∆lsym || ∆lt ≈ ∆lsym || throw(ArgumentError("∆lt = $∆lt must be greater than or equal to ∆lsym = $∆lsym."))
    (L = diff(gap)[1]) > 0 || throw(ArgumentError("gap[2] = $(gap[2]) must be strictly greater than gap[1] = $(gap[1])."))

    if issmooth(SVector(∆lsym, ∆lt), rt)
        filler = fill_constant(∆lsym, gap, ∆lt, rt, rmax)
    else
        ∆lmax = ∆lt
        ∆lmin = ∆lsym

        # Guess graded ∆l's.
        n = ceil(log(∆lmax/∆lmin)/log(rt))  # smallest integer n satisfying r = (∆lmax/∆lmin)^(1/n) ≤ rt
        r = (∆lmax/∆lmin)^(1/n)  # fastest geometric ratio (satisfying ≤ rt) for sequence to grow from ∆lmin to ∆lmax
        ∆l_array = ∆lmin * (r.^(1:n))

        L_graded = sum(∆l_array)  # ∆lmin * (r^1 + ... + r^n)
        if 2L_graded > L  # filler ≈ [∆l_array; reverse(∆l_array)], hence the factor 2
            # In this case, ∆l will remain smaller than ∆lt.  Just increase ∆l from ∆lsym as
            # fast as possible (i.e., at the geometric ratio close to the target ratio rt)
            # and decrease it before it gets too large, such that the sum of ∆l's doesn't
            # spill over the gap.

            # Below, n[1] is such that ceil(n[1]) is the smallest integer satisfying
            # ∆lmin * (rt^1 + ... + rt^n + ... + rt^1) ≥ L; note rt^n is added once.
            # Similarly, n[2] is such that ceil(n[2]) is the smallest integer satisfying
            # ∆lmin * (rt^1 + ... + rt^n + rt^n ... + rt^1) ≥ L; note rt^n is added twice.
            ns = VecFloat(undef, 2)
            # n[1] = fzero(n -> ∆lmin*rt * (rt^(n-1) + 2(rt^(n-1) - 1)/(rt-1)) - L, 1)
            # n[2] = fzero(n -> ∆lmin*rt * 2(rt^n - 1)/(rt-1) - L, 1)
            ns[1], isconverged = newtsol(
                1.,
                n -> ∆lmin * rt * (rt^(n-1) + 2(rt^(n-1)-1)/(rt-1)) - L,
                n -> ∆lmin * rt^n * log(rt) * (rt+1) / (rt-1)
            )
            @assert isconverged
            ns[2], isconverged = newtsol(
                1.,
                n -> 2∆lmin * rt * (rt^n-1)/(rt-1) - L,
                n -> 2∆lmin * rt^(n+1) * log(rt) / (rt-1)
            )
            @assert isconverged

            # Between Float ns[1] and ns[2], choose the one closer to their smallest
            # upper-bounding integers.
            _m, i = findmin(abs.(ceil.(ns) - ns))
            n = ceil(ns[i])

            # ns[1] > ns[2], and ns[2] > 0 because f(n) = ∆lmin*rt * (rt^(n-1) + 2(rt^(n-1) - 1)/(rt-1))
            # is an increasing function of n and f(0) = -∆lmin < 0, whereas ns[2] is the
            # solution to  f(n) = L > 0.  Therefore, ceil(ns[1]) ≥ ceil(ns[2]) ≥ 1.
            @assert n ≥ 1

            if i == 1
                # Below, ∆lmin * (r^1 + ... + r^n + ... + r^1) == L
                # r = fzero(s -> ∆lmin*s * (s^(n-1) + 2*(s^(n-1) - 1)/(s-1)) - L, rt)
                r, isconverged = newtsol(
                    rt,
                    s -> ∆lmin * s * (s^(n-1) + 2*(s^(n-1) - 1)/(s-1)) - L,
                    s -> ∆lmin * (s^(n-1) * (n*s^2 - 2s - n) + 2) / (s-1)^2
                )
                @assert isconverged
                ∆l_array = ∆lmin * (r.^[1:n; n-1:-1:1])
            else  # i == 2
                @assert i==2
                # Below, ∆lmin * (r^1 + ... + r^n + r^n + ... + r^1) == L
                # r = fzero(s -> ∆lmin*s * (s^n - 1)/(s-1) - L/2, rt)
                r, isconverged = newtsol(
                    rt,
                    s -> ∆lmin * s * (s^n-1)/(s-1) - L/2,
                    s -> ∆lmin * (s^n * (n*s - n - 1) + 1) / (s-1)^2
                )
                @assert isconverged
                ∆l_array = ∆lmin * (r.^[1:n; n:-1:1])
            end
            @assert r ≤ rt
            r ≥ 1/rt || throw(ArgumentError("gap = $gap is too small for ∆lsym = $∆lsym to grow geometrically to ∆lt = $∆lt with geometric ratio ≤ rt = $rt."))

            ∆l_filler = ∆l_array
        else  # 2L_graded ≤ L
            # Slightly underfill the gap with ∆lmax's plus the above generated graded ∆l's.
            n_∆lmax = floor(Int, (L - 2L_graded) / ∆lmax)  # use floor() to underfill
            ∆lmax_array = fill(∆lmax, n_∆lmax)  # [∆lmax, ..., ∆lmax]
            L_∆lmax = sum(∆lmax_array)

            # Slightly overfill the gap with ∆lmin's, graded ∆l's, and the above generated L_∆lmax.
            @assert 2L_graded + L_∆lmax ≤ L
            n_∆lmin = ceil(Int, (L - 2L_graded - L_∆lmax) / 2∆lmin)  # use ceil() to overfill
            ∆lmin_array = fill(∆lmin, n_∆lmin)  # [∆lmin, ..., ∆lmin]
            L_∆lmin = sum(∆lmin_array)

            # Reduce ∆l's slightly such that they fit into the space remaining after filling L with ∆lmax's and ∆lmin's.
            L_graded = (L - L_∆lmax - 2L_∆lmin) / 2
            # Below, ∆lmin * (rnew^1 + ... + rnew^n) == L_graded
            # rnew = fzero(s -> ∆lmin*s * (s^n - 1) / (s-1) - L_graded, r)
            rnew, isconverged = newtsol(
                r,
                s -> ∆lmin*s * (s^n - 1) / (s-1) - L_graded,
                s -> ∆lmin * (s^n * (n*s - n - 1) + 1) / (s-1)^2
            )
            @assert isconverged
            @assert rnew ≤ r
            r = rnew
            ∆l_array = ∆lmin * (r.^(1:n))
            ∆l_filler = [∆lmin_array; ∆l_array; ∆lmax_array; reverse(∆l_array); ∆lmin_array]
        end

        # Construct filler.
        # If we set filler = cumsum([gap[1]; ∆l_filler]), the last entry of filler
        # might not be precisely gap[2] due to round-off error, so do the following.
        filler_fwd = cumsum([gap[1]; ∆l_filler])  # last element ≈ gap[2]
        filler_bwd = cumsum([gap[2]; -reverse(∆l_filler)])  # last element ≈ gap[1]
        filler = [gap[1]; (@view(filler_fwd[2:end-1]) + reverse(@view(filler_bwd[2:end-1]))) / 2; gap[2]]
    end

    return filler
end

function fill_geometric(∆ln::Real, gap::AbsVecReal, ∆lt::Real, ∆lp::Real, rt::Real=R_TARGET_DEFAULT, rmax::Real=R_MAX_DEFAULT)
    # Generate a filler subgrid whose grid vertex separations increases
    # geometrically from both ends.

    # gap: length-2 row vector whose elements are the locations of the beginning and
    # ending grid vertices of the gap
    # ∆ln, ∆lp: ∆l of the subgrids before and after the gap
    # ∆lt: target ∆l in the gap
    # rt: target ratio of geometric sequence
    # rmax: maximum ratio of geometric sequence

    ((∆lt ≥ ∆ln || ∆lt ≈ ∆ln) && (∆lt ≥ ∆lp || ∆lt ≈ ∆lp)) ||
        throw(ArgumentError("∆lt = $∆lt must be greater than or equal to ∆ln = $∆ln and ∆lp = $∆lp."))
    (L = diff(gap)[1]) > 0 || throw(ArgumentError("gap[2] = $(gap[2]) must be strictly greater than gap[1] = $(gap[1])."))

    if ∆ln ≈ ∆lp
        filler = fill_geometric_sym(∆ln, gap, ∆lt, rt, rmax)
    else  # ∆ln < ∆lp or ∆ln > ∆lp
        ∆lmax = max(∆ln, ∆lp)
        ∆lmin = min(∆ln, ∆lp)

        # Guess graded ∆l's.
        n = ceil(log(∆lmax/∆lmin)/log(rt))  # smallest integer n satisfying r = (∆lmax/∆lmin)^(1/n) ≤ rt
        r = (∆lmax/∆lmin)^(1/n)  # fastest geometric ratio (satisfying ≤ rt) for sequence to grow frow ∆lmin to ∆lmax
        ∆l_array = ∆lmin * (r.^(1:n))

        # Slightly underfill the gap with ∆lmax and the above generated graded ∆l's.
        L_graded = sum(∆l_array)  # ∆lmin * (r^1 + ... + r^n)
        L_graded ≤ L || throw(ArgumentError("gap = $gap is too small for min(∆ln,∆lp) = $(min(∆ln,∆lp)) to grow geometrically to max(∆ln,∆lp) = $(max(∆ln,∆lp)) with geometric ratio ≤ rt = $rt."))
        if ∆ln < ∆lp
            filler_n = cumsum([gap[1]; ∆l_array])
            gap_sym = [filler_n[end]; gap[2]]
            filler_sym = fill_geometric_sym(∆lp, gap_sym, ∆lt, rt, rmax)
            filler = [filler_n[1:end-1]; filler_sym]
        else
            @assert ∆ln > ∆lp
            filler_p = cumsum([gap[2]; -∆l_array])
            filler_p = reverse(filler_p)
            gap_sym = [gap[1]; filler_p[1]]
            filler_sym = fill_geometric_sym(∆ln, gap_sym, ∆lt, rt, rmax)
            filler = [filler_sym; filler_p[2:end]]
        end
    end

    return filler
end

"""
    comp_lprim1d

Complete a 1D primal grid by filling the gaps between parts of the grid.

`comp_lprim1d(sublprim)` attempts to generate a full 1D primal
grid from a partially generated 1D primal grid (called "subgrids").

`sublprim` is an array whose each entry is an array of Float.  An example of `sublprim`
is `[[0., 0.5], [1.], [10., 11., 12.], [1.5], [20.5, 21.], [2.], [30., 31.]]`.
The entries at odd-numbered indices are in the form of `[l₁, ..., lₙ]`, and they
represent a subgrid (a complete grid between `l₁` and `lₙ`).  Between two adjacent
subgrids, at even-numbered indices of `sublprim`, come entries in the form of `[∆lt]`,
and they represent the target grid cell size between two subgrids.  `comp_lprim1d`
tries to complete grid generation, by filling the spaces between subgrids with grid
cells of size `∆lt`.

In the above example, `[0., 0.5]`, `[10., 11., 12.]`, `[20.5, 21.]`, `[30., 31.]`
are subgrids; `[1.]` is the target `∆l` between `0.5` and `10.`; `1.5` is the target
`∆l` between `12.` and `20.5`; `2.` is the target `∆l` between `21.` and `30.`.
Then, `comp_lprim1d(sublprim)` generates a grid between `0.` and `31.`.  Errors
are generated when the function fails to generate a primal grid.

```jldoctest
julia> sublprim = [[0., 0.5], [1.], [10., 11., 12], [1.5], [20.5, 21.], [2.], [30., 31.]];
julia> lprim = comp_lprim1d(sublprim);
```
"""
function comp_lprim1d(sublprim::AbsVec{<:AbsVecReal}, rt::Real=R_TARGET_DEFAULT, rmax::Real=R_MAX_DEFAULT)
    # Check inputs.
    nentry = length(sublprim)
    if nentry == 1
        lprim = sublprim[1]
        return lprim
    end

    # Now nentry ≥ 2.
    isodd(nentry) || throw(ArgumentError("sublprim = $sublprim must have odd number of entries."))

    i = 1
    sublprim_curr = sublprim[i]  # [l₁, ..., lₙ]
    length(sublprim_curr) ≥ 2 || throw(ArgumentError("Entry #$(i) of sublprim, $(sublprim_curr), must be length-2 or longer."))
    issorted(sublprim_curr) || throw(ArgumentError("Entry #$(i) of sublprim, $(sublprim_curr), must be sorted in ascending order."))

    lprim = copy(sublprim_curr)
    for i = 3:2:nentry
        sublprim_next = sublprim[i]  # [l₁, ..., lₙ]
        ∆lt = sublprim[i-1]  # [∆l]

        length(sublprim_next) ≥ 2 || throw(ArgumentError("Entry #$(i) of sublprim, $(sublprim_next), must be length-2 or longer."))
        issorted(sublprim_next) || throw(ArgumentError("Entry #$(i) of sublprim, $(sublprim_next), must be sorted in ascending order."))
        length(∆lt) == 1 || throw(ArgumentError("Entry #$(i-1) of sublprim, $(∆lt), must be length-1."))
        sublprim_curr[end] < sublprim_next[1] || throw(ArgumentError("Last entry of subgrid $(sublprim_curr) must be strictly less than first entry of subgrid $(sublprim_next) of sublprim."))

        ∆ln = lprim[end] - lprim[end-1]  # == sublprim_curr[end] - sublprim_curr[end-1]; ∆l of current subgrid
        ∆lp = sublprim_next[2] - sublprim_next[1]  # ∆l of next subgrid
        gap = [lprim[end], sublprim_next[1]]  # gap range between neighboring subgrids
        try
            filler = fill_geometric(∆ln, gap, ∆lt[1], ∆lp, rt, rmax)

            # Below, note that the provided subgrids are preserved despite round-off
            # errors in `filler`.  This is important for assiging materials and
            # sources at intended locations.  For example, failure to maintain the
            # provided subgrids can result in warnings in Shape.generate_kernel().
            @assert filler[1]≈sublprim_curr[end] && filler[end]≈sublprim_next[1]
            append!(lprim, @view(filler[2:end-1]))
            append!(lprim, sublprim_next)

            sublprim_curr = sublprim_next
        catch err
            isa(err, ArgumentError) ?
                throw(ErrorException("Grid generation failed between subgrids $(sublprim_curr) and $(sublprim_next) with target ∆l = $(∆lt): " * err.msg)) :
                throw(err)
        end

        sublprim_prev = sublprim_curr
    end

    ∆lprim = diff(lprim)
    ind = findfirst_stiff_∆∆l(∆lprim, rmax)
    ind == 0 || throw(ErrorException("Grid generation failed: grid nodes $(lprim[ind:ind+2])
        are separated by $(∆lprim[ind:ind+1]), which do not vary smoothly."))

    return lprim
end


function gen_lprim1d(domain::OpenInterval,
                     domaintype::GridType,
                     intvprim::AbsVec{ClosedInterval},
                     intvdual::AbsVec{ClosedInterval},
                     lprim0::AbsVecReal,
                     ldual0::AbsVecReal,
                     rt::Real=R_TARGET_DEFAULT,
                     rmax::Real=R_MAX_DEFAULT)
    # intvprim = unique(intvprim)
    # lprim0 = unique(lprim0)
    # ldual0 = unique(ldual0)

    sublprim = gen_sublprim1d(domain, domaintype, intvprim, intvdual, lprim0, ldual0)
    lprim = comp_lprim1d(sublprim, rt, rmax)


    # Npml(w,Sign.n) = length(findall(lprim < lprim(1) + Lpml(w,Sign.n)));
    # Npml(w,Sign.p) = length(findall(lprim > lprim(end) - Lpml(w,Sign.p)));

    return lprim
end

function gen_lprim(domain::Object{N,<:Box{N}},
                   domaintype::NTuple{N,GridType},
                   sprim::AbsVec{<:Object{N}},
                   lprim0::NTuple{N,AbsVecReal}=ntuple(n->Float[],Val(N)),
                   ldual0::NTuple{N,AbsVecReal}=ntuple(n->Float[],Val(N)),
                   rt::Real=R_TARGET_DEFAULT,
                   rmax::Real=R_MAX_DEFAULT) where {N}
    lprim = [gen_lprim1d(OpenInterval(domain,w), domaintype[w], (s->ClosedInterval(s,w)).(sprim),
                         ClosedInterval[], lprim0[w], ldual0[w], rt, rmax) for w = 1:N]

    return ntuple(w->lprim[w], Val(N))
end

# function gen_lprim3d(domain::Interval3D,
#                      domaintype::Tuple3{GridType},
#                      intvprim::AbsVec{Interval3D},
#                      lprim0::Tuple3{AbsVecReal},
#                      ldual0::Tuple3{AbsVecReal},
#                      rt::Real=R_TARGET_DEFAULT,
#                      rmax::Real=R_MAX_DEFAULT)
#     xprim = genlprim1d(domain.comp[nX], domaintype[nX], (i->i.comp[nX]).(intvprim), lprim0[nX], ldual[nX])
#     yprim = genlprim1d(domain.comp[nY], domaintype[nY], (i->i.comp[nY]).(intvprim), lprim0[nY], ldual[nY])
#     zprim = genlprim1d(domain.comp[nZ], domaintype[nZ], (i->i.comp[nZ]).(intvprim), lprim0[nZ], ldual[nZ])
#
#     return (xprim, yprim, zprim)
# end

# function gen_lprim3d(domain::Interval3D, intervals::AbsArr{Interval3D}, withuniform::Bool = false)
#     lprims = cell(1, Axis.count);
#     Npml = NaN(Axis.count, Sign.count);
#     if withuniform
#         for w = Axis.elems
#             dl_intended = domain.dl_max(w);
#             Nw = round(domain.L(w) / dl_intended);
#             lprim = range(domain.bound(w,Sign.n), stop=domain.bound(w,Sign.p), length=Nw+1);
#             @assert length(lprim) >= 2
#             dl_generated = lprim(2) - lprim(1);
#             error_dl = (dl_generated - dl_intended) / dl_intended;
#             if abs(error_dl) > 1e-9
#                 warning('Maxwell:gridGen', ['grid vertex separation %s in generated uniform grid', ...
#                     'significantly differs from intended separation %s: error is %s percent.'], ...
#                     num2str(dl_generated), num2str(dl_intended), num2str(error_dl * 100));
#             end
#             Npml(w,Sign.n) = length(findall(lprim < lprim(1) + Lpml(w,Sign.n)));
#             Npml(w,Sign.p) = length(findall(lprim > lprim(end) - Lpml(w,Sign.p)));
#             lprims{w} = lprim;
#         end
#     else  # withuniform == false: use adaptive grid generation algorithm
#         # Below, all for loops over Axis.elems could be merged into a single for
#         # loop.  However, that is less efficient because, for example, it makes a
#         # shape retrieved three times instead of once from shapes.  When
#         # shapes is long, this turns out to be actually quite a large penalty.
#
#         intervals = cell(1, Axis.count);
#         for j = 1:length(shapes)
#             shape = shapes(j);
#             for w = Axis.elems
#                 inters_w = intervals{w};  # initially empty
#                 inter_curr = shape.interval(w);
#                 is_new = true;
#                 i = 1;
#                 while is_new && i <= length(inters_w)
#                     is_new = is_new && ~isequal(inters_w(i), inter_curr);  # isequal() compares contents of two objects
#                     i = i + 1;
#                 end
#
#                 # Keep only a new interval; many intervals can be the same, e.g., in a
#                 # photonic crystal.
#                 if is_new
#                     intervals{w} = [inters_w, inter_curr];
#                 end
#             end
#         end
#
#         lprim0 = cell(1, Axis.count);
#         ldual0 = cell(1, Axis.count);
#         for j = 1:length(srcs)
#             src = srcs(j);
#             for w = Axis.elems
#                 lprim0{w} = [lprim0{w}, src.l{w, GT.prim}];
#                 ldual0{w} = [ldual0{w}, src.l{w, GT.dual}];
#             end
#         end
#
#         for w = Axis.elems
#             try
#                 lprim_part = generate_lprim1d_part(domain.interval(w), Lpml(w,:), intervals{w}, lprim0{w}, ldual0{w});
#                 lprim = complete_lprim1d(lprim_part);
#             catch err1
#                 try
#                     lprim_part = generate_lprim1d_part(domain.interval(w), Lpml(w,:), intervals{w}, lprim0{w}, []);
#                     lprim = complete_lprim1d(lprim_part);
#                     ldual = mean([lprim(1:end-1); lprim(2:end)]);  # take average in column
#                 catch err2
#                     exception = MException('Maxwell:gridGen', '%s-axis grid generation failed.', char(w));
#                     throw(addCause(exception, err2));
#                 end
#
#                 if ~isempty(setdiff(ldual0{w}, ldual))
#                     exception = MException('Maxwell:gridGen', '%s-axis grid generation failed.', char(w));
#                     throw(addCause(exception, err1));
#                 end
#             end
#             Npml(w,Sign.n) = length(findall(lprim < lprim(1) + Lpml(w,Sign.n)));
#             Npml(w,Sign.p) = length(findall(lprim > lprim(end) - Lpml(w,Sign.p)));
#             lprims{w} = lprim;
#         end
#     end
#
#     return (lprims, Npml)
# end
