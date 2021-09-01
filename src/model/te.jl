# Alias of Model{K,Kₑ,Kₘ, K², K₊₁,K₊₂, AK₊₁,AK₊₂} for 2D TE Maxwell's equations.
# The type parameters AK₊₁ and AK₊₂ specify device-specific arrays types (e.g., CuArray) and
# are user-defined in the constructor.
const ModelTE{AK₊₁,AK₊₂} = Model{2,2,1, 4, 3,4, AK₊₁,AK₊₂}

# Convenience constructor
function ModelTE(grid::Grid; Atype::Type=Array)
    cmpₛ = SInt(1,2)  # shapes in xy-plane
    cmpₑ = SInt(1,2)  # E-field with x- and y-components
    cmpₘ = SInt(3)  # H-field with z-component

    return ModelTE{Atype{ComplexF,3},Atype{ComplexF,4}}(;grid, cmpₛ, cmpₑ, cmpₘ)
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
    gₑ = ft2gt.(EE,boundft)
    gₘ = ft2gt.(HH,boundft)

    εarr = mdl.εarr
    µarr = mdl.μarr

    oind2shp = mdl.oind2shp
    oind2εind = mdl.oind2εind
    oind2μind = mdl.oind2μind

    cmpₑ = mdl.cmpₑ
    cmpₘ = mdl.cmpₘ

    εind2ε = sub_pind2matprm(mdl.εind2ε, cmpₑ)
    μind2μ = sub_pind2matprm(mdl.μind2μ, cmpₘ)

    # Create temporary storages.
    εxx_oind2d = create_oind_array(N)
    εyy_oind2d = create_oind_array(N)
    εoo_oind2d = create_oind_array(N)

    μzz_oind2d = create_oind_array(N)

    # Assign the material parameters.
    assign_param!(εarr, (εyy_oind2d,εxx_oind2d), oind2shp, oind2εind, εind2ε, gₑ, τl, isbloch)  # diagonal entries of ε tensors
    assign_param!(εarr, tuple(μzz_oind2d), oind2shp, oind2εind, εind2ε, gₑ, τl, isbloch)  # off-diagonal entries of ε tensors
    assign_param!(μarr, tuple(εoo_oind2d), oind2shp, oind2μind, μind2μ, gₘ, τl, isbloch)  # μ tensors (rank-0, so scalars)

    # Smooth the material parameters.
    smooth_param!(εarr, (εxx_oind2d,εyy_oind2d), oind2shp, oind2εind, εind2ε, gₑ, l, lg, σ, ∆τ, mdl.ise˔shp)  # diagonal entries of ε tensors
    smooth_param!(εarr, tuple(εoo_oind2d), oind2shp, oind2εind, εind2ε, gₑ, l, lg, σ, ∆τ, mdl.ise˔shp)  # off-diagonal entries of ε tensors
    smooth_param!(μarr, tuple(μzz_oind2d), oind2shp, oind2μind, μind2μ, gₘ, l, lg, σ, ∆τ, mdl.ish˔shp)  # μ tensors (rank-0, so scalars)

    return nothing
end
