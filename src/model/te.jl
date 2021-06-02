# Type abbreviation of Model{K,Kₑ,Kₘ,K²,Kₑ²,Kₘ²,K₊₁,K₊₂} for 2D TE Maxwell's equations
const ModelTE = Model{2,2,1, 4,4,1, 3,4}

# Convenience constructor
function ModelTE(grid::Grid)
    cmpₛ = SVec(1,2)  # shapes in xy-plane
    cmpₑ = SVec(1,2)  # E-field with x- and y-components
    cmpₘ = SVec(3)  # H-field with z-component

    return ModelTE(grid, cmpₛ, cmpₑ, cmpₘ)
end

# Assign the material parameters on the grid and smooth them.
function calc_matparams!(mdl::ModelTE)
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
    εxx_oind2d = create_oind_array(N)
    εyy_oind2d = create_oind_array(N)
    εoo_oind2d = create_oind_array(N)

    μzz_oind2d = create_oind_array(N)

    # Assign the material parameters.
    assign_param!(εarr, (εyy_oind2d,εxx_oind2d), oind2shp, oind2εind, εind2ε, ge, τl, isbloch)  # diagonal entries of ε tensors
    assign_param!(εarr, tuple(μzz_oind2d), oind2shp, oind2εind, εind2ε, ge, τl, isbloch)  # off-diagonal entries of ε tensors
    assign_param!(μarr, tuple(εoo_oind2d), oind2shp, oind2μind, μind2μ, gh, τl, isbloch)  # μ tensors (rank-0, so scalars)

    # Smooth the material parameters.
    smooth_param!(εarr, (εxx_oind2d,εyy_oind2d), oind2shp, oind2εind, εind2ε, ge, l, lg, σ, ∆τ)  # diagonal entries of ε tensors
    smooth_param!(εarr, tuple(εoo_oind2d), oind2shp, oind2εind, εind2ε, ge, l, lg, σ, ∆τ)  # off-diagonal entries of ε tensors
    smooth_param!(μarr, tuple(μzz_oind2d), oind2shp, oind2μind, μind2μ, gh, l, lg, σ, ∆τ)  # μ tensors (rank-0, so scalars)

    return nothing
end
