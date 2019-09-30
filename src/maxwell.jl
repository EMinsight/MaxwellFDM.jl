# Here, I am assuming that the E-field grid is primal in all the x-, y-, z-directions.  This
# is wrong.  Whether it is primal or dual can change depending on the boundary condition;
# see the leading comments in grid.jl.
#
# I have made the same mistake in assignment.jl and smoothing.jl.  I need to fix these all
# together someday.

# The getters like get_Amatrix actually change the contents of Maxwell by creating the
# desired fields on the fly.  However, because they do not change the contents of Maxwell
# in the ways noticeable to the users, I leave out the exclamation mark at the end of the
# function names.
#
# (I think I had to implement this on-the-fly creation in getters because some fields are
# not ready from the beginning.  Instead, I could choose to implement initialize! function
# that needs to be called after all the setters are called.  This is still an option, but I
# think there were cases some getters need to be called before all the setters are called.)
export Maxwell
export set_unitlen!, set_bounds!, set_∆l!, set_isbloch!, set_kbloch!, set_Npml!, set_wvlen!,
    set_freq!, get_unit, get_osc, get_grid, set_background!, add_obj!, get_param3d,
    get_stretched_∆l, get_εmatrix, get_curle, get_curlm, get_curls, get_dblcurl, get_Amatrix,
    get_Mc, get_Ml, get_Mr, add_srce!, add_srcm!, get_bvector

# Add quantities, and construct various systems at the end at once?
# Create a domain from the domain size, and add it to the object list.

mutable struct Maxwell
    # Unit
    unitlen::Real

    # Oscillation
    λ₀::Number

    # Grid
    g::Grid{3}
    bounds::Tuple2{AbsVecReal}  # ([xmin,ymin,zmin], [xmax,ymax,zmax])
    ∆l::AbsVecReal  # [∆x,∆y,∆z]
    boundft::SVec3{FieldType}  # boundary field type
    isbloch::AbsVecBool  # [Bool,Bool,Bool]
    kbloch::AbsVecReal  # [Real,Real,Real]
    e⁻ⁱᵏᴸ::AbsVecComplex

    # PML
    Npml::Tuple2{AbsVecInteger}
    s∆l::Tuple2{Tuple3{Vector{CFloat}}}

    # Domain
    εdom::Real

    # Objects and materials
    ovec::AbsVec{Object3}
    paramset::Tuple2{AbsVec{SMat3Complex}}
    param3d::Tuple2{AbsArrComplex{5}}

    # Sources
    je3d::AbsArrComplex{4}
    jm3d::AbsArrComplex{4}

    # Matrices and right-hand-side vector
    Mε::SparseMatrixCSC{CFloat,Int}
    Ce::SparseMatrixCSC{CFloat,Int}
    Cm::SparseMatrixCSC{CFloat,Int}
    CC::SparseMatrixCSC{CFloat,Int}
    A::SparseMatrixCSC{CFloat,Int}
    Mc::SparseMatrixCSC{CFloat,Int}  # operator interpolating Ex, Ey, Ez at grid cell corners
    Ml::SparseMatrixCSC{CFloat,Int}  # operator interpolating voxel-corner Ex, Ey, Ez at Ez, Ex, Ey (i.e., left component) locations
    Mr::SparseMatrixCSC{CFloat,Int}  # operator interpolating voxel-corner Ex, Ey, Ez at Ey, Ez, Ex (i.e., right component) locations
    b::Vector{CFloat}

    function Maxwell()
        m = new()

        # Initialize some fields with default values.
        m.boundft = SVector(EE,EE,EE)
        m.isbloch = @SVector ones(Bool, 3)
        m.kbloch = @SVector zeros(3)
        m.∆l = @SVector ones(3)

        m.ovec = Object3[]
        m.paramset = (SMat3Complex[], SMat3Complex[])

        return m
    end
end

#= Setters for basic quantities =#
set_unitlen!(m::Maxwell, unitlen::Real) = (m.unitlen = unitlen; return nothing)
set_bounds!(m::Maxwell, bounds::Tuple2{AbsVecReal}) = (m.bounds = bounds; return nothing)
set_∆l!(m::Maxwell, ∆l::AbsVecReal) = (m.∆l = ∆l; return nothing)
set_boundft(m::Maxwell, boundft::AbsVec{FieldType}) = (m.∆boundft = boundft; return nothing)
set_isbloch!(m::Maxwell, isbloch::AbsVecBool) = (m.isbloch = isbloch; return nothing)
set_kbloch!(m::Maxwell, kbloch::AbsVecReal) = (m.kbloch = kbloch; return nothing)
set_Npml!(m::Maxwell, Npml::Tuple2{AbsVecInteger}) = (m.Npml = Npml; return nothing)
set_wvlen!(m::Maxwell, λ₀::Number) = (m.λ₀ = λ₀; return nothing)
set_freq!(m::Maxwell, ω₀::Number) = (m.λ₀ = 2π/ω₀; return nothing)

#= Getters for basic constructed objects =#
get_unit(m::Maxwell) = PhysUnit(m.unitlen)
get_osc(m::Maxwell) = Oscillation(m.λ₀, get_unit(m))

function get_grid(m::Maxwell)
    if ~isdefined(m, :g)
        L = m.bounds[nP] - m.bounds[nN]
        N = round.(Int, L ./ m.∆l)
        lprim = map((lmin,lmax,n)->collect(range(lmin, stop=lmax, length=n+1)), m.bounds[nN], m.bounds[nP], N)

        m.g = Grid(get_unit(m), (lprim...,), m.isbloch)
    end

    return m.g
end

function get_e⁻ⁱᵏᴸ(m::Maxwell)
    if ~isdefined(m, :e⁻ⁱᵏᴸ)
        if iszero(m.kbloch)
            m.e⁻ⁱᵏᴸ = @SVector(ones(3))
        else
            g = get_grid(m)
            L = g.L
            kbloch = SVec3(m.kbloch)

            m.e⁻ⁱᵏᴸ = exp.(-im .* kbloch .* L)
        end
    end

    return m.e⁻ⁱᵏᴸ
end

#= Setters for objects =#
# Set background materials.
set_background!(m::Maxwell, matname::String, ε::MatParam) = add_obj!(m, matname, ε, Box(get_grid(m).bounds))

# Below, I write two methods for add_obj: one for a tuple and the other for a vector.
# Because we can easily create a tuple from a vector and vice versa, we can implement only
# one and make the other simply a wrapper of that.  Because transforming a tuple to a vector
# is more efficient than the opposite, I implement add_obj for a vector.
add_obj!(m::Maxwell, matname::String, ε::MatParam, shapes::Shape...) = add_obj!(m, matname, ε, [shapes...])

function add_obj!(m::Maxwell, matname::String, ε::MatParam, shapes::AbsVec{<:Shape})
    mat = Material(matname, ε=ε)
    for s = shapes  # shapes is tuple
        obj = Object(s, mat)
        add!(m.ovec, m.paramset, obj)
    end

    return nothing
end

#= Getters for operators =#
function get_param3d(m::Maxwell)
    if ~isdefined(m, :param3d)
        g = get_grid(m)
        N = g.N

        # Initialize other fields that depend on the grid.
        ε3d = create_param3d(N)
        εobj3d = create_n3d(Object3, N)
        εind3d = create_n3d(ParamInd, N)
        εoind3d = create_n3d(ObjInd, N)

        assign_param!(ε3d, εobj3d, εind3d, εoind3d, m.ovec, g.ghosted.τl, g.isbloch, m.boundft)
        smooth_param!(ε3d, εobj3d, εind3d, εoind3d, g.l, g.ghosted.l, g.σ, g.ghosted.∆τ, m.boundft)

        μ3d = create_param3d(N)
        μobj3d = create_n3d(Object3, N)
        μind3d = create_n3d(ParamInd, N)
        μoind3d = create_n3d(ObjInd, N)

        assign_param!(μ3d, μobj3d, μind3d, μoind3d, m.ovec, g.ghosted.τl, g.isbloch, m.boundft)
        smooth_param!(μ3d, μobj3d, μind3d, μoind3d, g.l, g.ghosted.l, g.σ, g.ghosted.∆τ, m.boundft)

        m.param3d = (ε3d, μ3d)
    end

    return m.param3d
end

function get_stretched_∆l(m::Maxwell)
    if ~isdefined(m, :s∆l)
        g = get_grid(m)
        lpml, Lpml = get_pml_loc(g.l[nPR], g.bounds, m.Npml)

        ω = in_ω₀(get_osc(m))
        sfactor = gen_stretch_factor(ω, g.l, lpml, Lpml)

        s∆lprim = map((x,y)->x.*y, sfactor[nPR], g.∆l[nPR])
        s∆ldual = map((x,y)->x.*y, sfactor[nDL], g.∆l[nDL])

        m.s∆l = (s∆lprim, s∆ldual)
    end

    return m.s∆l
end

function get_εmatrix(m::Maxwell)
    if ~isdefined(m, :Mε)
        g = get_grid(m)
        s∆l = get_stretched_∆l(m)
        param3d = get_param3d(m)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m)

        m.Mε = param3d2mat(param3d[nPR], [PRIM,PRIM,PRIM], g.N, s∆l[nDL], s∆l[nPR], g.isbloch, e⁻ⁱᵏᴸ, reorder=true)
    end

    return m.Mε
end

function get_curle(m::Maxwell)
    if ~isdefined(m, :Ce)
        g = get_grid(m)
        s∆l = get_stretched_∆l(m)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m)

        m.Ce = create_curl(m.boundft.==EE, g.N, s∆l[nDL], g.isbloch, e⁻ⁱᵏᴸ, reorder=true)
    end

    return m.Ce
end

function get_curlm(m::Maxwell)
    if ~isdefined(m, :Cm)
        g = get_grid(m)
        s∆l = get_stretched_∆l(m)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m)

        m.Cm = create_curl(m.boundft.==HH, g.N, s∆l[nPR], g.isbloch, e⁻ⁱᵏᴸ, reorder=true)
    end

    return m.Cm
end

get_curls(m::Maxwell) = (get_curle(m), get_curlm(m))

function get_dblcurl(m::Maxwell)
    if ~isdefined(m, :CC)
        Ce, Cm = get_curls(m)

        # Cm and Ce are sparse, but if they have explicit zeros, Cm * I drops those zeros
        # from Cm, leading to the symbolic sparsity pattern that is not the transpose of the
        # symbolic sparsity pattern of Ce (which still has explicit zeros).  This makes the
        # symbolic sparsity pattern of CC nonsymmetric, for which UMFPACK uses a slower
        # LU factorization algorithm.
        #
        # One way to avoid this problem is to use a sparse diagonal matrix instead of I for
        # Tμ⁻¹.  Another way is to drop the explicit zeros in Ce and Cm such that there are
        # no explicit zeros to drop in Ce and Cm.  I chose the latter, because it leads to a
        # compact CC that is faster to factorize.

        Tμ⁻¹ = I
        # g = get_grid(m)
        # M = 3*prod(g.N)
        # Tμ⁻¹ = sparse(Diagonal(ones(M)))

        m.CC = Cm * Tμ⁻¹ * Ce
    end

    return m.CC
end

function get_Amatrix(m::Maxwell)
    if ~isdefined(m, :A)
        CC = get_dblcurl(m)
        Mε = get_εmatrix(m)
        ω = in_ω₀(get_osc(m))

        m.A = CC - ω^2 * Mε;
    end

    return m.A
end

# Create the operator interpolating Ex, Ey, Ez at grid cell corners.
function get_Mc(m::Maxwell)
    if ~isdefined(m, :Mc)
        g = get_grid(m)
        s∆l = get_stretched_∆l(m)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m)

        # Arguments of create_mean:
        # - isfwd = (m.boundft.==HH) because we are doing backward averaging in the
        # w-direction if the w-boundary field is the E-field boundary, in order to get the
        # field at the voxel corner locations from the E-field locations.
        # - We supply ∆l and ∆l′ because we want to use weighted arithmetic averaging for
        # this backward averaging.
        # - kdiag = 0 for putting Ex, Ey, Ez to voxel corners.
        m.Mc = create_mean(m.boundft.==HH, g.N, s∆l[nDL], s∆l[nPR], g.isbloch, e⁻ⁱᵏᴸ, kdiag=0, reorder=true)
    end

    return m.Mc
end

# Create operator interpolating the voxel-corner Ex, Ey, Ez at Ez, Ex, Ey (i.e., left
# component locations).
function get_Ml(m::Maxwell)
    if ~isdefined(m, :Ml)
        g = get_grid(m)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m)

        # Arguments of create_mean:
        # - isfwd = (m.boundft.==EE) because we are doing forward averaging in the
        # w-direction if the w-boundary field is the E-field boundary, in order to get the
        # field at the E-field locations from the voxel corner locations.
        # - We do not supply ∆l and ∆l′ because we want to use unweighted arithmetic
        # averaging for this forward averaging.
        # - kdiag = +1 for putting Ex, Ey, Ez to Ez, Ex, Ey locations.
        m.Ml = create_mean(m.boundft.==EE, g.N, g.isbloch, e⁻ⁱᵏᴸ, kdiag=1, reorder=true)
    end

    return m.Ml
end

# Create operator interpolating the voxel-corner Ex, Ey, Ez at Ey, Ez, Ex (i.e., right
# component) locations.
function get_Mr(m::Maxwell)
    if ~isdefined(m, :Mr)
        g = get_grid(m)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m)

        # Arguments of create_mean:
        # - isfwd = (m.boundft.==EE) because we are doing forward averaging in the
        # w-direction if the w-boundary field is the E-field boundary, in order to get the
        # field at the E-field locations from the voxel corner locations.
        # - We do not supply ∆l and ∆l′ because we want to use unweighted arithmetic
        # averaging for this forward averaging.
        # - kdiag = -1 for putting Ex, Ey, Ez to Ey, Ez, Ex locations.
        m.Mr = create_mean(m.boundft.==EE, g.N, g.isbloch, e⁻ⁱᵏᴸ, kdiag=-1, reorder=true)
    end

    return m.Mr
end

#= Setters and getters for sources =#
function add_srce!(m::Maxwell, src::Source)
    if ~isdefined(m, :je3d)
        g = get_grid(m)
        m.je3d = create_field3d(g.N)
    end

    add!(m.je3d, PRIM, g.bounds, g.l, g.∆l, g.isbloch, src)

    return nothing
end

function add_srcm!(m::Maxwell, src::Source)
    if ~isdefined(m, :jm3d)
        g = get_grid(m)
        m.jm3d = create_field3d(g.N)
    end

    add!(m.jm3d, DUAL, g.bounds, g.l, g.∆l, g.isbloch, src)

    return nothing
end

function get_bvector(m::Maxwell)
    if ~isdefined(m, :b)
        je = field3d2vec(m.je3d, reorder=true)
        jm = field3d2vec(m.jm3d, reorder=true)

        ω = in_ω₀(get_osc(m))
        Cm = get_curlm(m)
        m.b = -im * ω * je - Cm * jm
    end

    return m.b
end
