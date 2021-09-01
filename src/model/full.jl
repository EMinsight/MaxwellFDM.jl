# Alias of Model{K,Kₑ,Kₘ, K², K₊₁,K₊₂, AK₊₁,AK₊₂} for full 3D Maxwell's equations.
# The type parameters AK₊₁ and AK₊₂ specify device-specific arrays types (e.g., CuArray) and
# are user-defined in the constructor.
const ModelFull{AK₊₁,AK₊₂} = Model{3,3,3, 9, 4,5, AK₊₁,AK₊₂}

# Convenience constructor
function ModelFull(grid::Grid; Atype::Type=Array)
    cmpₛ = SInt(1,2,3)  # shapes in 3D space
    cmpₑ = SInt(1,2,3)  # E-field with x-, y-, and y-components
    cmpₘ = SInt(1,2,3)  # H-field with x-, y-, and y-components

    return ModelFull{Atype{ComplexF,4},Atype{ComplexF,5}}(;grid, cmpₛ, cmpₑ, cmpₘ)
end

# Assign the material parameters on the grid and smooth them.
function calc_matparams!(mdl::ModelFull)
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
    εxx_oind3d = create_oind_array(N)
    εyy_oind3d = create_oind_array(N)
    εzz_oind3d = create_oind_array(N)
    εoo_oind3d = create_oind_array(N)

    μxx_oind3d = create_oind_array(N)
    μyy_oind3d = create_oind_array(N)
    μzz_oind3d = create_oind_array(N)
    μoo_oind3d = create_oind_array(N)

    # Assign the material parameters.
    assign_param!(εarr, (μxx_oind3d,μyy_oind3d,μzz_oind3d), oind2shp, oind2εind, εind2ε, gₑ, τl, isbloch)  # diagonal entries of ε tensors
    assign_param!(εarr, tuple(μoo_oind3d), oind2shp, oind2εind, εind2ε, gₑ, τl, isbloch)  # off-diagonal entries of ε tensors
    assign_param!(μarr, (εxx_oind3d,εyy_oind3d,εzz_oind3d), oind2shp, oind2μind, μind2μ, gₘ, τl, isbloch)  # diagonal entries of μ tensors
    assign_param!(μarr, tuple(εoo_oind3d), oind2shp, oind2μind, μind2μ, gₘ, τl, isbloch)  # off-diagonal entries of μ tensors

    # Smooth the material parameters.
    smooth_param!(εarr, (εxx_oind3d,εyy_oind3d,εzz_oind3d), oind2shp, oind2εind, εind2ε, gₑ, l, lg, σ, ∆τ, mdl.ise˔shp)  # diagonal entries of ε tensors
    smooth_param!(εarr, tuple(εoo_oind3d), oind2shp, oind2εind, εind2ε, gₑ, l, lg, σ, ∆τ, mdl.ise˔shp)  # off-diagonal entries of ε tensors
    smooth_param!(μarr, (μxx_oind3d,μyy_oind3d,μzz_oind3d), oind2shp, oind2μind, μind2μ, gₘ, l, lg, σ, ∆τ, mdl.ish˔shp)  # diagonal entries of μ tensors
    smooth_param!(μarr, tuple(μoo_oind3d), oind2shp, oind2μind, μind2μ, gₘ, l, lg, σ, ∆τ, mdl.ish˔shp)  # off-diagonal entries of μ tensors

    return nothing
end
