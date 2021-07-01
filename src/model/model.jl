export set_ω!, set_boundft!, set_Npml!, set_kbloch!  # basic setter functions
export clear_objs!, clear_srcs!  # plural because clearing both electric and magnetic quantities
export add_srce!, add_srcm!  # add_obj! is imported from MaxwellBase
export create_paramops, create_curls, create_srcs, create_linsys, create_A, create_b, e2h, h2e
export create_Mcs, create_Mls, create_Mrs

# Do not export Model; quality it with the package name MaxwellFDFD, because I would have
# similar types in other packages such as MaxwellSALT and MaxwellGuide.
mutable struct Model{K,Kₑ,Kₘ,K₊₁,K₊₂,
                     AK₊₁<:AbsArrComplexF{K₊₁},AK₊₂<:AbsArrComplexF{K₊₂},
                     K²,Kₑ²,Kₘ²}
    # Frequency
    ω::Number  # can be complex

    # Dimensions
    grid::Grid{K}
    cmpₛ::SInt{K}  # Cartesian components of shape dimension
    cmpₑ::SInt{Kₑ}  # Cartesian components of E-field
    cmpₘ::SInt{Kₘ}  # Cartesian components of H-field

    # Boundary properties
    boundft::SVec{K,FieldType}
    Npml::Tuple2{SInt{K}}
    kbloch::SFloat{K}  # [kx_bloch, ky_bloch, kz_bloch]

    # Material parameter arrays
    εarr::AK₊₂
    μarr::AK₊₂

    # Current density arrays
    jₑarr::AK₊₁
    jₘarr::AK₊₁

    # Storage for assignment and smoothing of material parameters
    oind2shp::Vector{Shape{K,K²}}
    oind2εind::Vector{ParamInd}
    oind2μind::Vector{ParamInd}
    εind2ε::Vector{SSComplexF{Kₑ,Kₑ²}}
    μind2μ::Vector{SSComplexF{Kₘ,Kₘ²}}

    function Model{K,Kₑ,Kₘ,K₊₁,K₊₂,AK₊₁,AK₊₂,K²,Kₑ²,Kₘ²}(
        grid::Grid{K}, cmpₛ::SInt{K}, cmpₑ::SInt{Kₑ}, cmpₘ::SInt{Kₘ}
        ) where {K,Kₑ,Kₘ,K₊₁,K₊₂,AK₊₁<:AbsArrComplexF{K₊₁},AK₊₂<:AbsArrComplexF{K₊₂},K²,Kₑ²,Kₘ²}
        ω = 0.0

        εarr = create_param_array(grid.N, ncmp=Kₑ)  # filled with zeros
        µarr = create_param_array(grid.N, ncmp=Kₘ)  # filled with zeros

        jₑarr = create_field_array(grid.N, ncmp=Kₑ)  # filled with zeros
        jₘarr = create_field_array(grid.N, ncmp=Kₘ)  # filled with zeros

        oind2shp = Shape{K,K^2}[]
        oind2εind = ParamInd[]
        oind2μind = ParamInd[]
        εind2ε = SSComplexF{Kₑ,Kₑ^2}[]
        μind2μ = SSComplexF{Kₘ,Kₘ^2}[]

        boundft = SVec(ntuple(k->EE, Val(K)))
        kbloch = SVec(ntuple(k->0.0, Val(K)))
        Npml = (SVec(ntuple(k->0, Val(K))), SVec(ntuple(k->0, Val(K))))

        return new(ω, grid, cmpₛ, cmpₑ, cmpₘ, boundft, Npml, kbloch, εarr, μarr, jₑarr, jₘarr,
                   oind2shp, oind2εind, oind2μind, εind2ε, μind2μ)
    end
end

# Basic setters
set_ω!(mdl::Model{K}, ω::Number) where {K} = (mdl.ω = ω; nothing)
set_boundft!(mdl::Model{K}, boundft::AbsVec{FieldType}) where {K} = (mdl.boundft = SVec{K}(boundft); nothing)
set_Npml!(mdl::Model{K}, Npml::Tuple2{AbsVecInteger}) where {K} = (mdl.Npml = SVec{K}.(Npml); nothing)
set_kbloch!(mdl::Model{K}, kbloch::AbsVecReal) where {K} = (mdl.kbloch = SVec{K}(kbloch); nothing)

create_e⁻ⁱᵏᴸ(mdl::Model) = exp.(-im .* mdl.kbloch .* mdl.grid.L)

# Main functions
function clear_objs!(mdl::Model{K,Kₑ,Kₘ}) where {K,Kₑ,Kₘ}
    mdl.εarr .= 0
    mdl.μarr .= 0

    empty!(mdl.oind2shp)
    empty!(mdl.oind2εind)
    empty!(mdl.oind2μind)
    empty!(mdl.εind2ε)
    empty!(mdl.μind2μ)

    return nothing
end

MaxwellBase.add_obj!(mdl::Model{K}, matname::String, shapes::AbsVec{<:Shape{K}};
                     ε::MatParam=1.0, μ::MatParam=1.0) where {K} =
    add_obj!(mdl, matname, ε=ε, μ=μ, tuple(shapes...))

function MaxwellBase.add_obj!(mdl::Model{K,Kₑ,Kₘ}, matname::String, shapes::Shape{K}...;
                              ε::MatParam=1.0, μ::MatParam=1.0) where {K,Kₑ,Kₘ}
    mat = Material{Kₑ,Kₘ}(matname, ε=ε, μ=μ)
    for shp = shapes  # shapes is tuple
        obj = Object(shp, mat)
        add_obj!(mdl.oind2shp, (mdl.oind2εind,mdl.oind2μind), (mdl.εind2ε,mdl.μind2μ), obj)
    end

    return nothing
end

function create_stretched_∆ls(mdl::Model)
    ω = mdl.ω
    grid = mdl.grid
    Npml = mdl.Npml
    s∆l = create_stretched_∆l(ω, grid, Npml)
    s∆l⁻¹ = invert_∆l(s∆l)

    boundft = mdl.boundft
    gₑ = ft2gt.(EE, boundft)
    gₘ = ft2gt.(HH, boundft)

    s∆lₑ = t_ind(s∆l, gₑ)  # stretched ∆l's centered at E-field plane locations
    s∆lₘ = t_ind(s∆l, gₘ)  # stretched ∆l's centered at H-field plane locations
    s∆lₑ⁻¹ = t_ind(s∆l⁻¹, gₑ)
    s∆lₘ⁻¹ = t_ind(s∆l⁻¹, gₘ)

    return s∆lₑ, s∆lₘ, s∆lₑ⁻¹, s∆lₘ⁻¹
end

function create_paramops(mdl::Model; order_cmpfirst::Bool=true)
    s∆lₑ, s∆lₘ, s∆lₑ⁻¹, s∆lₘ⁻¹ = create_stretched_∆ls(mdl)
    calc_matparams!(mdl)

    boundft = mdl.boundft
    gₑ = ft2gt.(EE, boundft)
    gₘ = ft2gt.(HH, boundft)

    isbloch = mdl.grid.isbloch
    e⁻ⁱᵏᴸ = create_e⁻ⁱᵏᴸ(mdl)
    Mε = create_paramop(mdl.εarr, gₑ, s∆lₘ, s∆lₑ⁻¹, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)
    Mμ = create_paramop(mdl.μarr, gₘ, s∆lₑ, s∆lₘ⁻¹, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)

    return Mε, Mμ
end

function create_curls(mdl::Model; order_cmpfirst::Bool=true)
    _, _, s∆lₑ⁻¹, s∆lₘ⁻¹ = create_stretched_∆ls(mdl)

    boundft = mdl.boundft
    cmpₛ, cmpₑ, cmpₘ = mdl.cmpₛ, mdl.cmpₑ, mdl.cmpₘ

    isbloch = mdl.grid.isbloch
    e⁻ⁱᵏᴸ = create_e⁻ⁱᵏᴸ(mdl)
    Cₑ = create_curl(boundft.==EE, s∆lₘ⁻¹, isbloch, e⁻ⁱᵏᴸ; cmp_shp=cmpₛ, cmp_out=cmpₘ, cmp_in=cmpₑ, order_cmpfirst)
    Cₘ = create_curl(boundft.==HH, s∆lₑ⁻¹, isbloch, e⁻ⁱᵏᴸ; cmp_shp=cmpₛ, cmp_out=cmpₑ, cmp_in=cmpₘ, order_cmpfirst)

    return Cₑ, Cₘ
end

function clear_srcs!(mdl::Model)
    mdl.jₑarr .= 0
    mdl.jₘarr .= 0

    return nothing
end

function add_srce!(mdl::Model{K,Kₑ}, src::Source{K,Kₑ}) where {K,Kₑ}
    gₑ = ft2gt.(EE, mdl.boundft)
    grid = mdl.grid
    jₑarr = mdl.jₑarr
    add_src!(jₑarr, gₑ, grid.bounds, grid.l, grid.∆l, grid.isbloch, src)

    return nothing
end

function add_srcm!(mdl::Model{K,Kₑ,Kₘ}, src::Source{K,Kₘ}) where {K,Kₑ,Kₘ}
    gₘ = ft2gt.(HH, mdl.boundft)
    grid = mdl.grid
    jₘarr = mdl.jₘarr
    add_src!(jₘarr, gₘ, grid.bounds, grid.l, grid.∆l, grid.isbloch, src)

    return nothing
end

function create_srcs(mdl::Model; order_cmpfirst::Bool=true)
    jₑ = field_arr2vec(mdl.jₑarr; order_cmpfirst)
    jₘ = field_arr2vec(mdl.jₘarr; order_cmpfirst)

    return jₑ, jₘ  # note jₑ and jₘ do not share memory with mdl.jₑarr and mdl.jₘarr
end

create_linsys(ft::FieldType, ω::Number, M::Tuple2{AbsMatNumber}, C::Tuple2{AbsMatNumber}, j::Tuple2{AbsVecNumber}) =
    create_linsys(ft, ω, M..., C..., j...)

function create_linsys(ft::FieldType, ω::Number,
                       Mε::AbsMatNumber, Mμ::AbsMatNumber,
                       Cₑ::AbsMatNumber, Cₘ::AbsMatNumber,
                       jₑ::AbsVecNumber, jₘ::AbsVecNumber)
    A = create_A(ft, ω, Mε, Mμ, Cₑ, Cₘ)
    b = create_b(ft, ω, Mε, Mμ, Cₑ, Cₘ, jₑ, jₘ)

    return A, b
end

create_A(ft::FieldType, ω::Number, M::Tuple2{AbsMatNumber}, C::Tuple2{AbsMatNumber}) =
    create_A(ft, ω, M..., C...)

function create_A(ft::FieldType,
                  ω::Number,
                  Mε::AbsMatNumber, Mμ::AbsMatNumber,
                  Cₑ::AbsMatNumber, Cₘ::AbsMatNumber)
    # Consider a potential approach to save memory:
    #     Y = copy(Cₑ)
    #     ldiv!(some_factorization(Mμ), Y)
    #     A = copy(Y)
    #     mul!(A, Cₘ, Y)
    # This still uses two allocations, so it is no better than A = Cₘ * (Mμ \ Cₑ)
    if ft == EE  # A = Cₘ (Mμ⁻¹ Cₑ) - ω² Mε
        A = Cₘ * (Mμ \ Cₑ)
        iszero(ω) || (A -= ω^2 * Mε)
    elseif ft == HH  # A = Cₑ (Mε⁻¹ Cₘ) - ω² Mμ
        A = Cₑ * (Mε \ Cₘ)
        iszero(ω) || (A -= ω^2 * Mμ)
    else
        @error "ft = $ft is unsupported."
    end

    return A
end

create_b(ft::FieldType, ω::Number, M::Tuple2{AbsMatNumber}, C::Tuple2{AbsMatNumber}, j::Tuple2{AbsVecNumber}) =
    create_b(ft, ω, M..., C..., j...; order_cmpfirst)

function create_b(ft::FieldType, ω::Number,
                  Mε::AbsMatNumber, Mμ::AbsMatNumber,
                  Cₑ::AbsMatNumber, Cₘ::AbsMatNumber,
                  jₑ::AbsVecNumber, jₘ::AbsVecNumber)
    # Consider a potential approach to save memory:
    #     y = copy(jₘ)
    #     ldiv!(some_factorization(Mμ), y)
    #     b = copy(y)
    #     mul!(b, Cₘ, y)
    #     b .*= -1
    # This still uses two allocations, so it is no better than b = -Cₘ * (Mμ \ jₘ)
    if ft == EE  # b = -i ω jₑ - Cₘ (Mμ⁻¹ jₘ)
        b = Cₘ * (Mμ \ jₘ)
        b .*= -1
        iszero(ω) || (b .-= (im * ω) .* jₑ)
    elseif ft == HH  # b = -i ω jₘ + Cₑ (Mε⁻¹ jₑ)
        b = Cₑ * (Mε \ jₑ)
        iszero(ω) || (b .-= (im * ω) .* jₘ)
    else
        @error "ft = $ft is unsupported."
    end

    return b
end

e2h(e::AbsVecNumber, ω::Number, M::Tuple2{AbsMatNumber}, C::Tuple2{AbsMatNumber}, j::Tuple2{AbsVecNumber}) =
    e2h(e, ω, M[2], C[1], j[2])
e2h(e::AbsVecNumber, ω::Number, Mμ::AbsMatNumber, Cₑ::AbsMatNumber, jₘ::AbsVecNumber) =
    (im/ω) .* (Mμ \ (Cₑ*e + jₘ))

h2e(h::AbsVecNumber, ω::Number, M::Tuple2{AbsMatNumber}, C::Tuple2{AbsMatNumber}, j::Tuple2{AbsVecNumber}) =
    h2e(h, ω, M[1], C[2], j[1])
h2e(h::AbsVecNumber, ω::Number, Mε::AbsMatNumber, Cₘ::AbsMatNumber, jₑ::AbsVecNumber) =
    (-im/ω) .* (Mε \ (Cₘ*h - jₑ))

# Create the operator interpolating solution fields at grid cell corners.
function create_Mcs(mdl::Model; order_cmpfirst::Bool=true)
    s∆lₑ, s∆lₘ, s∆lₑ⁻¹, s∆lₘ⁻¹ = create_stretched_∆ls(mdl)

    # Arguments of create_mean:
    # - isfwd is set to average the solution fields to get the fields at the voxel corners.
    # (For boundft[w] = ft, solution field is backward-averaged along w.)
    # - ∆l and ∆l′⁻¹ are supplied to use weighted arithmetic averaging in calculating
    # the field at the voxel corners.
    boundft = mdl.boundft
    isbloch = mdl.grid.isbloch
    e⁻ⁱᵏᴸ = create_e⁻ⁱᵏᴸ(mdl)
    Mcₑ = create_mean(boundft.!=EE, s∆lₘ, s∆lₑ⁻¹, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)
    Mcₘ = create_mean(boundft.!=HH, s∆lₑ, s∆lₘ⁻¹, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)

    return Mcₑ, Mcₘ
end

# Create operator interpolating the voxel-corner solution fields Fx, Fy, Fz at Fz, Fx, Fy
# (i.e., left component locations).
#
# Note that the corner field averaged along the w-direction is put at the Fw-location.  For
# example, the corner Fx is averaged along the z-direction to be placed at the Fz-location.
function create_Mls(mdl::Model; order_cmpfirst::Bool=true)
    N = mdl.grid.N
    Pl = create_πcmp(N, SVec(3,1,2); order_cmpfirst)  # move Fx to Fz-location, Fy to Fx-location, and Fz to Fy-location

    # Arguments of create_mean:
    # - isfwd is set to average the fields at the voxel corners to get the fields at the
    # solution field locations (but with polarization not meant for that locations).
    # (For boundft[w] = ft, voxel-corner field is forward-averaged along w.)
    # - ∆l and ∆l′⁻¹ are not supplied; use unweighted arithmetic averaging.
    boundft = mdl.boundft
    isbloch = mdl.grid.isbloch
    e⁻ⁱᵏᴸ = create_e⁻ⁱᵏᴸ(mdl)
    Mₑ = create_mean(boundft.==EE, N, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)
    Mₘ = create_mean(boundft.==HH, N, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)

    Mlₑ = Mₑ * Pl
    Mlₘ = Mₘ * Pl

    return Mlₑ, Mlₘ
end

# Create operator interpolating the voxel-corner solution fields Fx, Fy, Fz at Fy, Fz, Fx
# (i.e., right component) locations.
#
# Note that the corner field averaged along the w-direction is put at the Fw-location.  For
# example, the corner Fx is averaged along the z-direction to be placed at the Fz-location.
function create_Mrs(mdl::Model; order_cmpfirst::Bool=true)
    N = mdl.grid.N
    Pr = create_πcmp(N, SVec(2,3,1); order_cmpfirst)  # move Fx to Fy-location, Fy to Fz-location, and Fz to Fx-location

    # Arguments of create_mean:
    # - isfwd is set to average the fields at the voxel corners to get the fields at the
    # solution field locations (but with polarization not meant for that locations).
    # - ∆l and ∆l′⁻¹ are not supplied; use unweighted arithmetic averaging.
    boundft = mdl.boundft
    isbloch = mdl.grid.isbloch
    e⁻ⁱᵏᴸ = create_e⁻ⁱᵏᴸ(mdl)
    Mₑ = create_mean(boundft.==EE, N, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)
    Mₘ = create_mean(boundft.==HH, N, isbloch, e⁻ⁱᵏᴸ; order_cmpfirst)

    Mrₑ = Mₑ * Pr
    Mrₘ = Mₘ * Pr

    return Mrₑ, Mrₘ
end

include("full.jl")
include("te.jl")
include("tm.jl")
include("tem.jl")
