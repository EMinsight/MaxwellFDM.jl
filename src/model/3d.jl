# Here, I am assuming that the E-field grid is primal in all the x-, y-, z-directions.  This
# is wrong.  Whether it is primal or dual can change depending on the boundary condition;
# see the leading comments in grid.jl.
#
# I have made the same mistake in assignment.jl and smoothing.jl.  I need to fix these all
# together someday.

# The getters like get_Amatrix actually change the contents of Model3D by creating the
# desired fields on the fly.  However, because they do not change the contents of Model3D
# in the ways noticeable to the users, I leave out the exclamation mark at the end of the
# function names.
#
# (I think I had to implement this on-the-fly creation in getters because some fields are
# not ready from the beginning.  Instead, I could choose to implement initialize! function
# that needs to be called after all the setters are called.  This is still an option, but I
# think there were cases some getters need to be called before all the setters are called.)
export set_unitlen!, set_bounds!, set_∆l!, set_boundft!, set_isbloch!, set_kbloch!, set_Npml!,
    set_wvlen!, set_freq!, get_unit, get_osc, get_grid, set_background!, add_obj!, get_param3d,
    get_stretched_∆l, get_εmatrix, get_μmatrix, get_curle, get_curlm, get_curls, get_dblcurl,
    get_Amatrix, get_Mc, get_Ml, get_Mr, add_srce!, add_srcm!, get_bvector
# export Model3D  # do not export; instantiate with qualification as MaxwellWave.Model3D(...)

# Add quantities, and construct various systems at the end at once?
# Create a domain from the domain size, and add it to the object list.

mutable struct Model3D
    # Unit
    unitlen::Float

    # Oscillation
    λ₀::CFloat

    # Grid
    g3::Grid{3}
    bounds::Tuple2{SFloat{3}}  # ([xmin,ymin,zmin], [xmax,ymax,zmax])
    ∆l::SFloat{3}  # [∆x,∆y,∆z]; assume uniform grid
    boundft::SVec{3,FieldType}  # [boundft_x, boundft_y, boundft_z]
    isbloch::SBool{3}  # [isbloch_x, isbloch_y, isbloch_z]
    kbloch::SFloat{3}  # [kx_bloch, ky_bloch, kz_bloch]
    e⁻ⁱᵏᴸ::SNumber{3}  # [exp(-i kx Lx), exp(-i ky Ly), exp(-i kz Lz)

    # PML
    Npml::Tuple2{SInt{3}}  # ([Npml_x⁻,Npml_y⁻,Npml_z⁻], [Npml_x⁺,Npml_y⁺,Npml_z⁺])
    s∆l::Tuple23{VecComplex}  # ((sx ∆xprim, sy ∆yprim, sz ∆zprim), (sx ∆xdual, sy ∆ydual, sz ∆zdual))
    s∆l⁻¹::Tuple23{VecComplex}  # ((sx⁻¹ ∆xprim⁻¹, sy⁻¹ ∆yprim⁻¹, sz⁻¹ ∆zprim⁻¹), (sx⁻¹ ∆xdual⁻¹, sy⁻¹ ∆ydual⁻¹, sz⁻¹ ∆zdual⁻¹))

    # Fields for objects assignment and smoothing
    oind2shp::Vector{Shape3}
    oind2εind::Vector{ParamInd}
    oind2μind::Vector{ParamInd}
    εind2ε::Vector{SSComplex3}
    μind2μ::Vector{SSComplex3}

    param3d::Tuple2{ArrComplex{5}}  # (ε3d, μ3d)

    # Sources
    je3d::ArrComplex{4}
    jm3d::ArrComplex{4}

    # Matrices and right-hand-side vector
    ftsol::FieldType  # EE or HH: field type of solution variable
    order_cmpfirst::Bool
    Mε::SparseMatrixCSC{CFloat,Int}
    Mμ::SparseMatrixCSC{CFloat,Int}
    Ce::SparseMatrixCSC{CFloat,Int}
    Cm::SparseMatrixCSC{CFloat,Int}
    CC::SparseMatrixCSC{CFloat,Int}
    A::SparseMatrixCSC{CFloat,Int}
    Mc::SparseMatrixCSC{CFloat,Int}  # operator interpolating Ex, Ey, Ez at grid cell corners
    Ml::SparseMatrixCSC{CFloat,Int}  # operator interpolating voxel-corner Ex, Ey, Ez at Ez, Ex, Ey (i.e., left component) locations
    Mr::SparseMatrixCSC{CFloat,Int}  # operator interpolating voxel-corner Ex, Ey, Ez at Ey, Ez, Ex (i.e., right component) locations
    b::Vector{CFloat}

    function Model3D()
        m3 = new()

        # Initialize some fields with default values.
        m3.boundft = SVec(EE,EE,EE)
        m3.isbloch = SVec(true,true,true)
        m3.kbloch = SVec(0.0, 0.0, 0.0)
        m3.e⁻ⁱᵏᴸ = SVec(1.0, 1.0, 1.0)

        m3.oind2shp = Shape3[]
        m3.oind2εind = ParamInd[]
        m3.oind2μind = ParamInd[]
        m3.εind2ε = SSComplex3[]
        m3.μind2μ = SSComplex3[]

        m3.ftsol = EE
        m3.order_cmpfirst = true

        return m3
    end
end

#= Setters for basic quantities =#
set_unitlen!(m3::Model3D, unitlen::Real) = (m3.unitlen = unitlen; return nothing)
set_bounds!(m3::Model3D, bounds::Tuple2{AbsVecReal}) = (m3.bounds = SFloat{3}.(bounds); return nothing)
set_∆l!(m3::Model3D, ∆l::AbsVecReal) = (m3.∆l = SFloat{3}(∆l); return nothing)
set_boundft!(m3::Model3D, boundft::AbsVec{FieldType}) = (m3.boundft = SVec{3,FieldType}(boundft); return nothing)
set_isbloch!(m3::Model3D, isbloch::AbsVecBool) = (m3.isbloch = SBool{3}(isbloch); return nothing)
set_kbloch!(m3::Model3D, kbloch::AbsVecReal) = (m3.kbloch = SFloat{3}(kbloch); return nothing)
set_Npml!(m3::Model3D, Npml::Tuple2{AbsVecInteger}) = (m3.Npml = SInt{3}.(Npml); return nothing)
set_wvlen!(m3::Model3D, λ₀::Number) = (m3.λ₀ = λ₀; return nothing)
set_freq!(m3::Model3D, ω₀::Number) = (m3.λ₀ = 2π/ω₀; return nothing)
set_order_cmpfirst!(m3::Model3D, order_cmpfirst::Bool) = (m3.order_cmpfirst = order_cmpfirst; return nothing)
set_ftsol!(m3::Model3D, ftsol::AbsVec{FieldType}) = (m3.ftsol = Vector(ftsol); return nothing)

#= Getters for basic constructed objects =#
get_unit(m3::Model3D) = PhysUnit(m3.unitlen)
get_osc(m3::Model3D) = Oscillation(m3.λ₀, get_unit(m3))

function get_grid(m3::Model3D)
    if ~isdefined(m3, :g3)
        L = m3.bounds[nP] - m3.bounds[nN]
        N = round.(Int, L ./ m3.∆l)
        lprim = map((lmin,lmax,n)->collect(range(lmin, stop=lmax, length=n+1)), m3.bounds[nN], m3.bounds[nP], N)

        m3.g3 = Grid((lprim...,), m3.isbloch)
    end

    return m3.g3
end

function get_e⁻ⁱᵏᴸ(m3::Model3D)
    if ~isdefined(m3, :e⁻ⁱᵏᴸ)
        if iszero(m3.kbloch)
            m3.e⁻ⁱᵏᴸ = @SVector(ones(3))
        else
            g3 = get_grid(m3)
            L = g3.L
            kbloch = SVec{3}(m3.kbloch)

            m3.e⁻ⁱᵏᴸ = exp.(-im .* kbloch .* L)
        end
    end

    return m3.e⁻ⁱᵏᴸ
end

#= Setters for objects =#
# Set background materials.
set_background!(m3::Model3D, matname::String; ε::MatParam=1.0, μ::MatParam=1.0) = add_obj!(m3, matname, ε=ε, μ=μ, Box(get_grid(m3).bounds))


# Below, I write two methods for add_obj: one for a tuple and the other for a vector.
# Because we can easily create a tuple from a vector and vice versa, we can implement only
# one and make the other simply a wrapper of that.  Because transforming a tuple to a vector
# is more efficient than the opposite, I implement add_obj for a vector.
#
# For the reason of qualification by MaxwellBase, see Agenda > MaxwellWave > Daily log on
# Jan/21/2021.
MaxwellBase.add_obj!(m3::Model3D, matname::String, shapes::AbsVec{<:Shape3}; ε::MatParam=1.0, μ::MatParam=1.0) =
    add_obj!(m3, matname, ε=ε, μ=μ, tuple(shapes...))

function MaxwellBase.add_obj!(m3::Model3D, matname::String, shapes::Shape3...; ε::MatParam=1.0, μ::MatParam=1.0)
    mat = Material{3,3}(matname, ε=ε, μ=μ)
    for s = shapes  # shapes is tuple
        obj = Object(s, mat)
        add_obj!(m3.oind2shp, (m3.oind2εind,m3.oind2μind), (m3.εind2ε,m3.μind2μ), obj)
    end

    return nothing
end

#= Getters for operators =#
function get_param3d(m3::Model3D)
    if ~isdefined(m3, :ε3d) || ~isdefined(m3, :μ3d)
        g3 = get_grid(m3)
        N = g3.N
        boundft = m3.boundft

        oind2shp = m3.oind2shp
        oind2εind = m3.oind2εind
        oind2μind = m3.oind2μind
        εind2ε = m3.εind2ε
        μind2μ = m3.μind2μ

        ε3d = create_param_array(N)
        εxx_oind3d = create_oind_array(N)
        εyy_oind3d = create_oind_array(N)
        εzz_oind3d = create_oind_array(N)
        εoo_oind3d = create_oind_array(N)

        μ3d = create_param_array(N)
        μxx_oind3d = create_oind_array(N)
        μyy_oind3d = create_oind_array(N)
        μzz_oind3d = create_oind_array(N)
        μoo_oind3d = create_oind_array(N)

        assign_param!(ε3d, (μxx_oind3d,μyy_oind3d,μzz_oind3d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g3.ghosted.τl, g3.isbloch)
        assign_param!(ε3d, tuple(μoo_oind3d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g3.ghosted.τl, g3.isbloch)
        assign_param!(μ3d, (εxx_oind3d,εyy_oind3d,εzz_oind3d), ft2gt.(HH,boundft), oind2shp, oind2μind, μind2μ, g3.ghosted.τl, g3.isbloch)
        assign_param!(μ3d, tuple(εoo_oind3d), ft2gt.(HH,boundft), oind2shp, oind2μind, μind2μ, g3.ghosted.τl, g3.isbloch)

        smooth_param!(ε3d, (εxx_oind3d,εyy_oind3d,εzz_oind3d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)
        smooth_param!(ε3d, tuple(εoo_oind3d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)
        smooth_param!(μ3d, (μxx_oind3d,μyy_oind3d,μzz_oind3d), oind2shp, oind2μind, μind2μ, ft2gt.(HH,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)
        smooth_param!(μ3d, tuple(μoo_oind3d), oind2shp, oind2μind, μind2μ, ft2gt.(HH,boundft), g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)

        m3.param3d = (ε3d, μ3d)
    end

    return m3.param3d
end

function get_stretched_∆l(m3::Model3D)
    if ~isdefined(m3, :s∆l)
        g3 = get_grid(m3)
        ω = in_ω₀(get_osc(m3))

        m3.s∆l = create_stretched_∆l(ω, g3, m3.Npml)
    end

    return m3.s∆l
end

function get_stretched_∆l⁻¹(m3::Model3D)
    if ~isdefined(m3, :s∆l⁻¹)
        m3.s∆l⁻¹ = invert_∆l(get_stretched_∆l(m3))
    end

    return m3.s∆l⁻¹
end

function get_εmatrix(m3::Model3D)
    if ~isdefined(m3, :Mε)
        g3 = get_grid(m3)
        s∆l = get_stretched_∆l(m3)
        s∆l⁻¹ = get_stretched_∆l⁻¹(m3)
        param3d = get_param3d(m3)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m3)

        ge = ft2gt.(EE, m3.boundft)  # ge[w]: grid type of locations of E-field orthogonal to w-direction
        gh = ft2gt.(HH, m3.boundft)  # gh[w]: grid type of locations of H-field orthogonal to w-direction
        ∆le = t_ind(s∆l, gh)  # Ew is multiplied with ∆w centered at location of H-field orthogonal to w-direction
        ∆lh⁻¹ = t_ind(s∆l⁻¹, ge)  # Ew is divided by ∆w′ centered at E-field voxel corner

        ε3d = m3.param3d[nE]
        m3.Mε = create_paramop(ε3d, ge, g3.N, ∆le, ∆lh⁻¹, g3.isbloch, e⁻ⁱᵏᴸ, order_cmpfirst=m3.order_cmpfirst)
    end

    return m3.Mε
end

function get_μmatrix(m3::Model3D)
    if ~isdefined(m3, :Mμ)
        g3 = get_grid(m3)
        s∆l = get_stretched_∆l(m3)
        s∆l⁻¹ = get_stretched_∆l⁻¹(m3)
        param3d = get_param3d(m3)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m3)

        gh = ft2gt.(HH, m3.boundft)  # gh[w]: grid type of locations of H-field orthogonal to w-direction
        ge = ft2gt.(EE, m3.boundft)  # ge[w]: grid type of locations of E-field orthogonal to w-direction
        ∆lh = t_ind(s∆l, ge)  # Hw is multiplied with ∆w centered at location of E-field orthogonal to w-direction
        ∆le⁻¹ = t_ind(s∆l⁻¹, gh)  # Hw is divided by ∆w′ centered at H-field voxel corner

        μ3d = m3.param3d[nH]
        m3.Mμ = create_paramop(μ3d, gh, g3.N, ∆lh, ∆le⁻¹, g3.isbloch, e⁻ⁱᵏᴸ, order_cmpfirst=m3.order_cmpfirst)
    end

    return m3.Mμ
end

function get_curle(m3::Model3D)
    if ~isdefined(m3, :Ce)
        g3 = get_grid(m3)
        s∆l⁻¹ = get_stretched_∆l⁻¹(m3)

        gh = ft2gt.(HH, m3.boundft)  # HH, not EE, because ∆w for Ce is centered at location of H-field orthogonal to w-direction
        ∆le⁻¹ = t_ind(s∆l⁻¹, gh)

        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m3)
        isfwd = m3.boundft.==EE
        m3.Ce = create_curl(isfwd, g3.N, ∆le⁻¹, g3.isbloch, e⁻ⁱᵏᴸ, order_cmpfirst=m3.order_cmpfirst)
    end

    return m3.Ce
end

function get_curlm(m3::Model3D)
    if ~isdefined(m3, :Cm)
        g3 = get_grid(m3)
        s∆l⁻¹ = get_stretched_∆l⁻¹(m3)

        ge = ft2gt.(EE, m3.boundft)  # EE, not HH, because ∆w for Cm is centered at location of E-field orthogonal to w-direction
        ∆lh⁻¹ = t_ind(s∆l⁻¹, ge)

        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m3)
        isfwd = m3.boundft.==HH
        m3.Cm = create_curl(isfwd, g3.N, ∆lh⁻¹, g3.isbloch, e⁻ⁱᵏᴸ, order_cmpfirst=m3.order_cmpfirst)
    end

    return m3.Cm
end

get_curls(m3::Model3D) = (get_curle(m3), get_curlm(m3))

function get_dblcurl(m3::Model3D)
    if ~isdefined(m3, :CC)
        Ce, Cm = get_curls(m3)

        # Notes about old code using Mμ⁻¹ = I
        #
        # Cm and Ce are sparse, but if they have explicit zeros, Cm * I drops those zeros
        # from Cm, leading to the symbolic sparsity pattern that is not the transpose of the
        # symbolic sparsity pattern of Ce (which still has explicit zeros).  This makes the
        # symbolic sparsity pattern of m3.CC nonsymmetric, for which UMFPACK uses a slower
        # LU factorization algorithm.
        #
        # One way to avoid this problem is to use, instead of I, a sparse diagonal matrix
        # for Tμ⁻¹ as
        #
        #     g3 = get_grid(m3)
        #     M = 3*prod(g3.N)
        #     Tμ⁻¹ = sparse(Diagonal(ones(M)))
        #
        # Another way is to drop the explicit zeros in Ce and Cm such that there are
        # no explicit zeros to drop in Ce and Cm.  I chose the latter, because it leads to a
        # compact CC that is faster to factorize.  This is implemented in create_curl().

        if m3.ftsol == EE
            Mμ = get_μmatrix(m3)
            m3.CC = Cm * (Mμ \ Ce)
        else  # m3.ftsol = HH
            Mε = get_εmatrix(m3)
            m3.CC = Ce * (Mε \ Cm)
        end
    end

    return m3.CC
end

function get_Amatrix(m3::Model3D)
    if ~isdefined(m3, :A)
        if m3.ftsol == EE
            CC = get_dblcurl(m3)
            Mε = get_εmatrix(m3)
            ω = in_ω₀(get_osc(m3))

            m3.A = CC - ω^2 * Mε;
        else # m3.ftsol = HH
            CC = get_dblcurl(m3)
            Mμ = get_μmatrix(m3)
            ω = in_ω₀(get_osc(m3))

            m3.A = CC - ω^2 * Mμ;
        end
    end

    return m3.A
end

# Create the operator interpolating solution fields at grid cell corners.
function get_Mc(m3::Model3D)
    if ~isdefined(m3, :Mc)
        g3 = get_grid(m3)
        s∆l = get_stretched_∆l(m3)
        s∆l⁻¹ = get_stretched_∆l⁻¹(m3)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m3)

        ft = m3.ftsol[1]  # field type of solution field
        ft′ = alter(ft)  # field type of complementary field
        gt = ft2gt.(ft, m3.boundft)  # voxel corners are at locations of solution field, whose field type is ft
        gt′ = ft2gt.(ft′, m3.boundft)  # voxel edges are at locations of field complementary to solution field

        isfwd = m3.boundft.==ft′  # for boundft[w] = ft, solution field is backward-averaged along w
        ∆l = t_ind(s∆l, gt′)  # multiply voxel edges
        ∆l′⁻¹ = t_ind(s∆l, gt)  # divide by complementary edges

        # Arguments of create_mean:
        # - isfwd is set to average the solution fields to get the fields at the voxel corners.
        # - ∆l and ∆l′⁻¹ are supplied to use weighted arithmetic averaging in calculating
        # the field at the voxel corners.
        m3.Mc = create_mean(isfwd, g3.N, ∆l, ∆l′⁻¹, g3.isbloch, e⁻ⁱᵏᴸ, order_cmpfirst=m3.order_cmpfirst)
    end

    return m3.Mc
end

# Create operator interpolating the voxel-corner solution fields Fx, Fy, Fz at Fz, Fx, Fy
# (i.e., left component locations).
#
# Note that the corner field averaged along the w-direction is put at the Fw-location.  For
# example, the corner Fx is averaged along the z-direction to be placed at the Fz-location.
# This means the m̂_w operator generates the w-component in the output, and therefore we pass
# outpermutem̂[w] = w to create_mean(), i.e., the keyword argument outpermutem̂ of create_mean()
# retains the default value.  On the other hand, inpermutem̂ must change from the default
# value; see the function body.
function get_Ml(m3::Model3D)
    if ~isdefined(m3, :Ml)
        g3 = get_grid(m3)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m3)

        ft = m3.ftsol[1]  # field type of solution field
        isfwd = m3.boundft.==ft  # for boundft[w] = ft, voxel-corner field is forward-averaged along w

        # Arguments of create_mean:
        # - isfwd is set to average the fields at the voxel corners to get the fields at the
        # solution field locations (but with polarization not meant for that locations).
        # - ∆l and ∆l′⁻¹ are not supplied; use unweighted arithmetic averaging.
        # - inpermutem̂ = [2,3,1] because m̂x, m̂y, m̂z are applied to voxel-corner Fy, Fz, Fx,
        # to put the averaged fields at the Fx-, Fy-, Fz-locations, respectively.
        m3.Ml = create_mean(isfwd, g3.N, g3.isbloch, e⁻ⁱᵏᴸ, inpermutem̂=SVec(2,3,1), order_cmpfirst=m3.order_cmpfirst)
    end

    return m3.Ml
end

# Create operator interpolating the voxel-corner solution fields Fx, Fy, Fz at Fy, Fz, Fx
# (i.e., right component) locations.
#
# Note that the corner field averaged along the w-direction is put at the Fw-location.  For
# example, the corner Fx is averaged along the z-direction to be placed at the Fz-location.
# This means the m̂_w operator generates the w-component in the output, and therefore we pass
# outpermutem̂[w] = w to create_mean(), i.e., the keyword argument outpermutem̂ of create_mean()
# retains the default value.  On the other hand, inpermutem̂ must change from the default
# value; see the function body.
function get_Mr(m3::Model3D)
    if ~isdefined(m3, :Mr)
        g3 = get_grid(m3)
        e⁻ⁱᵏᴸ = get_e⁻ⁱᵏᴸ(m3)

        ft = m3.ftsol[1]  # field type of solution field
        isfwd = m3.boundft.==ft  # for boundft[w] = ft, voxel-corner field is forward-averaged along w

        # Arguments of create_mean:
        # - isfwd is set to average the fields at the voxel corners to get the fields at the
        # solution field locations (but with polarization not meant for that locations).
        # - ∆l and ∆l′⁻¹ are not supplied; use unweighted arithmetic averaging.
        # - inpermutem̂ = [3,1,2] because m̂x, m̂y, m̂z are applied to voxel-corner Fz, Fx, Fy,
        # to put the averaged fields at the Fx-, Fy-, Fz-locations, respectively.
        m3.Mr = create_mean(isfwd, g3.N, g3.isbloch, e⁻ⁱᵏᴸ, inpermutem̂=SVec(3,1,2), order_cmpfirst=m3.order_cmpfirst)
    end

    return m3.Mr
end

#= Setters and getters for sources =#
function get_je3d(m3::Model3D)
    if ~isdefined(m3, :je3d)
        g3 = get_grid(m3)
        m3.je3d = create_field_array(g3.N)
    end

    return m3.je3d
end

function get_jm3d(m3::Model3D)
    if ~isdefined(m3, :jm3d)
        g3 = get_grid(m3)
        m3.jm3d = create_field_array(g3.N)
    end

    return m3.jm3d
end

function add_srce!(m3::Model3D, src::Source{3,3})
    je3d = get_je3d(m3)
    g3 = get_grid(m3)
    add_src!(je3d, ft2gt.(EE,m3.boundft), g3.bounds, g3.l, g3.∆l, g3.isbloch, src)

    return nothing
end

function add_srcm!(m3::Model3D, src::Source{3,3})
    jm3d = get_jm3d(m3)
    g3 = get_grid(m3)
    add_src!(jm3d, ft2gt.(HH,m3.boundft), g3.bounds, g3.l, g3.∆l, g3.isbloch, src)

    return nothing
end

function get_bvector(m3::Model3D)
    if ~isdefined(m3, :b)
        # To-do: this is where e⁻ⁱᵏᴸ needs to be applied.
        je3d = get_je3d(m3)
        jm3d = get_jm3d(m3)

        je = field_arr2vec(m3.je3d, order_cmpfirst=m3.order_cmpfirst)
        jm = field_arr2vec(m3.jm3d, order_cmpfirst=m3.order_cmpfirst)

        ω = in_ω₀(get_osc(m3))
        if m3.ftsol == EE
            Cm = get_curlm(m3)
            Mμ = get_μmatrix(m3)
            m3.b = -im * ω * je - Cm * (Mμ\jm)
        else  # m3.ftsol = HH
            Ce = get_curlm(m3)
            Mε = get_εmatrix(m3)
            m3.b = -im * ω * jm + Ce * (Mε\je)
        end

    end

    return m3.b
end
