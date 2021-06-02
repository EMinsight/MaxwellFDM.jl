# Type abbreviation of Model{K,Kₑ,Kₘ,K²,Kₑ²,Kₘ²,K₊₁,K₊₂} for 1D TEM Maxwell's equations
const ModelTEM = Model{1,1,1, 1,1,1, 2,3}

# Convenience constructor
function ModelTEM(grid::Grid)
    cmpₛ = SVec(3)  # shapes in z-axis
    cmpₑ = SVec(1)  # E-field with x-component
    cmpₘ = SVec(2)  # H-field with y-component

    return ModelTEM(grid, cmpₛ, cmpₑ, cmpₘ)
end

# Assign the material parameters on the grid and smooth them.
function calc_matparams!(mdl::ModelTEM)
    # Take necessary fields of mdl.
    grid = mdl.grid

    isbloch = grid.isbloch
    N = grid.N
    l = grid.l
    σ = grid.σ

    lg = grid.ghosted.l
    τl = grid.ghosted.τl
    ∆τ = grid.ghosted.∆τ

    boundft = mdl.boundft
    ge = ft2gt.(EE,boundft)
    gh = ft2gt.(HH,boundft)

    εarr = mdl.εarr
    µarr = mdl.μarr

    oind2shp = mdl.oind2shp
    oind2εind = mdl.oind2εind
    oind2μind = mdl.oind2μind
    εind2ε = mdl.εind2ε
    μind2μ = mdl.μind2μ

    # Create temporary storages.
    εxx_oind1d = create_oind_array(N)
    μyy_oind1d = create_oind_array(N)

    # Assign the material parameters.
    assign_param!(εarr, tuple(μyy_oind1d), oind2shp, oind2εind, εind2ε, ge, τl, isbloch)  # ε tensors (rand-0, so scalars)
    assign_param!(μarr, tuple(εxx_oind1d), oind2shp, oind2μind, μind2μ, gh, τl, isbloch)  # μ tensors (rand-0, so scalars)

    # Smooth the material parameters.
    smooth_param!(εarr, tuple(εxx_oind1d), oind2shp, oind2εind, εind2ε, ge, l, lg, σ, ∆τ)  # ε tensors (rand-0, so scalars)
    smooth_param!(μarr, tuple(μyy_oind1d), oind2shp, oind2μind, μind2μ, gh, l, lg, σ, ∆τ)  # μ tensors (rand-0, so scalars)

    return nothing
end
