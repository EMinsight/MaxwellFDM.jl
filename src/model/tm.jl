# Abbreviation of Model{K,Kₑ,Kₘ,K₊₁,K₊₂,AK₊₁,AK₊₂,K²,Kₑ²,Kₘ²} for 2D TM Maxwell's equations
const ModelTM{AK₊₁,AK₊₂} = Model{2,1,2, 3,4, AK₊₁,AK₊₂, 4,1,4}

# Convenience constructor
function ModelTM(grid::Grid; Atype::Type=Array)
    cmpₛ = SVec(1,2)  # shapes in xy-plane
    cmpₑ = SVec(3)  # E-field with z-component
    cmpₘ = SVec(1,2)  # H-field with x- and y-components

    return ModelTM{Atype{ComplexF,3},Atype{ComplexF,4}}(grid, cmpₛ, cmpₑ, cmpₘ)
end

# Assign the material parameters on the grid and smooth them.
function calc_matparams!(mdl::ModelTM)
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
    μxx_oind2d = create_oind_array(N)
    μyy_oind2d = create_oind_array(N)
    μoo_oind2d = create_oind_array(N)

    εzz_oind2d = create_oind_array(N)

    # Assign the material parameters.
    assign_param!(μarr, (μyy_oind2d,μxx_oind2d), oind2shp, oind2μind, μind2μ, gh, τl, isbloch)  # diagonal entries of μ tensors
    assign_param!(μarr, tuple(εzz_oind2d), oind2shp, oind2μind, μind2μ, gh, τl, isbloch)  # off-diagonal entries of μ tensors
    assign_param!(εarr, tuple(μoo_oind2d), oind2shp, oind2εind, εind2ε, ge, τl, isbloch)  # ε tensors (rank-0, so scalars)

    # Smooth the material parameters.
    smooth_param!(μarr, (μxx_oind2d,μyy_oind2d), oind2shp, oind2μind, μind2μ, gh, l, lg, σ, ∆τ)  # diagonal entries of μ tensors
    smooth_param!(μarr, tuple(μoo_oind2d), oind2shp, oind2μind, μind2μ, gh, l, lg, σ, ∆τ)  # off-diagonal entries of μ tensors
    smooth_param!(εarr, tuple(εzz_oind2d), oind2shp, oind2εind, εind2ε, ge, l, lg, σ, ∆τ)  # ε tensors (rank-0, so scalars)

    return nothing
end
