export Maxwell
export set_unitlen!, set_bounds!, set_∆l!, set_isbloch!, set_Npml!, set_wvlen!, set_freq!,
    get_unit, get_osc, get_grid, set_bgmat!, add_obj!, get_param3d, get_stretched_∆l,
    get_εmatrix, get_curle, get_curlm, get_curls, get_Amatrix, add_srce!, add_srcm!, get_bvector

# Add quantities, and construct various systems at the end at once?
# Create a domain from the domain size, and add it to the object list.

mutable struct Maxwell
    # Unit
    unitlen::Real

    # Oscillation
    λ₀::Real

    # Grid
    g::Grid{3}
    bounds::Tuple2{AbsVecReal}  # ([xmin,ymin,zmin], [xmax,ymax,zmax])
    ∆l::AbsVecReal  # [∆x,∆y,∆z]
    isbloch::AbsVec{Bool}  # [::Bool,::Bool,::Bool]
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
    A::SparseMatrixCSC{CFloat,Int}
    b::Vector{CFloat}

    function Maxwell()
        m = new()

        # Initialize some fields with default values.
        m.isbloch = @SVector ones(Bool, 3)
        m.e⁻ⁱᵏᴸ = @SVector ones(CFloat, 3)
        m.∆l = @SVector ones(3)

        m.ovec = Object3[]
        m.paramset = (SMat3Complex[], SMat3Complex[])

        return m
    end
end

#= Setters for basic quantities =#
set_unitlen!(m::Maxwell, unitlen::Real) = (m.unitlen = unitlen)
set_bounds!(m::Maxwell, bounds::Tuple2{AbsVecReal}) = (m.bounds = bounds)
set_∆l!(m::Maxwell, ∆l::AbsVecReal) = (m.∆l = ∆l)
set_isbloch!(m::Maxwell, isbloch::AbsVecBool) = (m.isbloch = isbloch)
set_Npml!(m::Maxwell, Npml::Tuple2{AbsVecInteger}) = (m.Npml = Npml)
set_wvlen!(m::Maxwell, λ₀::Real) = (m.λ₀ = λ₀)
set_freq!(m::Maxwell, ω₀::Real) = (m.λ₀ = 2π/ω₀)

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

#= Setters for objects =#
# Set background materials.
set_bgmat!(m::Maxwell, matname::String, ε::MatParam) = add_obj!(m, matname, ε, Box(get_grid(m).bounds))

# Below, I write two methods for add_obj: one for a tuple and the other for a vector.
# Because we can easily create a tuple from a vector and vice versa, we can implement only
# one and make the other simply a wrapper of that.  Because transforming a tuple to a vector
# is more efficient than the opposite, I implement add_obj for a vector.
add_obj!(m::Maxwell, matname::String, ε::MatParam, shapes::Shape...) = add_obj!(m, matname, ε, [shapes...])

function add_obj!(m::Maxwell, matname::String, ε::MatParam, shapes::AbsVec{<:Shape})
    mat = EncodedMaterial(PRIM, Material(matname, ε=ε))
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
        param3d = create_param3d(N)
        obj3d = create_n3d(Object3, N)
        pind3d = create_n3d(ParamInd, N)
        oind3d = create_n3d(ObjInd, N)

        assign_param!(param3d, obj3d, pind3d, oind3d, m.ovec, g.ghosted.τl, g.isbloch)
        smooth_param!(param3d, obj3d, pind3d, oind3d, g.l, g.ghosted.l, g.σ, g.ghosted.∆τ)

        m.param3d = param3d
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
        m.Mε = param3d2mat(param3d[nPR], [PRIM,PRIM,PRIM], g.N, s∆l[nDL], s∆l[nPR], g.isbloch, m.e⁻ⁱᵏᴸ, reorder=true)
    end

    return m.Mε
end

function get_curle(m::Maxwell)
    if ~isdefined(m, :Ce)
        g = get_grid(m)
        s∆l = get_stretched_∆l(m)

        m.Ce = create_curl([true,true,true], g.N, s∆l[nDL], g.isbloch, m.e⁻ⁱᵏᴸ, reorder=true)
    end

    return m.Ce
end

function get_curlm(m::Maxwell)
    if ~isdefined(m, :Cm)
        g = get_grid(m)
        s∆l = get_stretched_∆l(m)

        m.Cm = create_curl([false,false,false], g.N, s∆l[nPR], g.isbloch, m.e⁻ⁱᵏᴸ, reorder=true)
    end

    return m.Cm
end

get_curls(m::Maxwell) = (get_curle(m), get_curlm(m))

function get_Amatrix(m::Maxwell)
    if ~isdefined(m, :A)
        Mε = get_εmatrix(m)
        Tμ⁻¹ = I

        Ce, Cm = get_curls(m)

        ω = in_ω₀(get_osc(m))

        m.A =  Cm * Tμ⁻¹ * Ce - ω^2 * Mε;
    end

    return m.A
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
