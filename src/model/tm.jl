# Alias of Model{K,Kₑ,Kₘ, K²,Kₑ²,Kₘ², K₊₁,K₊₂, AK₊₁,AK₊₂} for 2D TM Maxwell's equations.
# The type parameters AK₊₁ and AK₊₂ specify device-specific arrays types (e.g., CuArray) and
# are user-defined in the constructor.
const ModelTM{AK₊₁,AK₊₂} = Model{2,1,2, 4,1,4, 3,4, AK₊₁,AK₊₂}

# Convenience constructor
function ModelTM(grid::Grid; Atype::Type=Array)
    cmpₛ = SInt(1,2)  # shapes in xy-plane
    cmpₑ = SInt(3)  # E-field with z-component
    cmpₘ = SInt(1,2)  # H-field with x- and y-components

    return ModelTM{Atype{ComplexF,3},Atype{ComplexF,4}}(;grid, cmpₛ, cmpₑ, cmpₘ)
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
    gₑ = ft2gt.(EE,boundft)
    gₘ = ft2gt.(HH,boundft)

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
    assign_param!(μarr, (μyy_oind2d,μxx_oind2d), oind2shp, oind2μind, μind2μ, gₘ, τl, isbloch)  # diagonal entries of μ tensors
    assign_param!(μarr, tuple(εzz_oind2d), oind2shp, oind2μind, μind2μ, gₘ, τl, isbloch)  # off-diagonal entries of μ tensors
    assign_param!(εarr, tuple(μoo_oind2d), oind2shp, oind2εind, εind2ε, gₑ, τl, isbloch)  # ε tensors (rank-0, so scalars)

    # Smooth the material parameters.
    smooth_param!(μarr, (μxx_oind2d,μyy_oind2d), oind2shp, oind2μind, μind2μ, gₘ, l, lg, σ, ∆τ, mdl.ish˔shp)  # diagonal entries of μ tensors
    smooth_param!(μarr, tuple(μoo_oind2d), oind2shp, oind2μind, μind2μ, gₘ, l, lg, σ, ∆τ, mdl.ish˔shp)  # off-diagonal entries of μ tensors
    smooth_param!(εarr, tuple(εzz_oind2d), oind2shp, oind2εind, εind2ε, gₑ, l, lg, σ, ∆τ, mdl.ise˔shp)  # ε tensors (rank-0, so scalars)

    return nothing
end
