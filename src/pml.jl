# About PML setup
#
# - Let's define the computational domain as the region where we perform computation, and
# the physical domain where actual physical phenomena are simulated.  The physical domain is
# the computational domain outside PML, because PML is an artificial region that is
# introduced to absorb physical waves.
#
# Currently, we specify the size of the computational domain first, and then specify PML
# such that it "eats in" the computational domain to leave the physical domain.  This means
# you always need to define the computational domain larger than the physical domain by the
# PML thickness.
#
# This seems inconvenient.  It may seem more convenient to specify the physical domain first,
# specify the PML thickness, and then construct the computational domain by attaching PML
# to the physical domain.  That way, we would be able to turn on and off the absorbing
# boundary without changing the size of the physical domain of interest.  Also, increasing
# the PML thickness would be done simply by increasing the PML thickness, wherease with the
# current method it requires increasing the computational domain as well to keep the size of
# the physical domain the same.
#
# However, there is a reason for the current method.  One of the most popular use cases of
# PML is the simulation of an infinite waveguide, where we need to define the waveguide
# inside PML.  Therefore, we cannot "hide" PML from the users: we have to let the users put
# objects inside PML.  In that case, it is better to let the users specify the size of the
# region including PML, because that way the users would know how deep they should place the
# waveguide.
#
# You way think that you hide PML from the users by automatically constructing the material
# parameters inside PML, by taking the material parameters on the boundary of the physical
# domain and extending it homogeneously into PML.  However, such a homogeneous extension is
# not always what we want.  For example, to construct PML for a phothonic crystal waveguide,
# we have to specify the photonic crystal structure inside PML.  In that case the material
# composition inside PML is not simply the homogeneous extension of their values at the PML
# boundary
#
# - Another question that could be raised about the current method is the way to set the
# PML thickness.  Currently we specify the number of grid cells inside PML and multiply the
# cell size ∆l to figure out the PML thickness, rather than specifying the PML thickness and
# divide it by ∆l to figure out the number of grid cells inside PML.
#
# The latter look like a better method, because the reflectance of PML is a function of the
# PML loss parameter and thickness.  Therefore, for a given PML loss parameter, the
# reflectance does not change with the spatial resolution of the grid if the PML thickness
# is maintained.  Conversely, if we fix the number of grid cells inside PML, the PML
# thickness changes with the spacial resolution, so the PML loss parameter needs to change
# as well to achieve the same PML reflectance.
#
# However, in my experience, when using iterative methods to solve Maxwell's equations, the
# number of grid cells inside PML affects the convergence of iterative methods significantly,
# regardless of the spatial resolution.  For example, when I doubled the spatial resolution
# while keeping the PML thickness and therefore doubled the number of grid cells inside PML,
# the iterative methods no longer converged.  On the other hand, even if I doubled the
# spatial resolution, if I keep the number of grid cells inside PML to halve the PML
# thickness, I still got convergence.
#
# Of course, I have to use a stronger PML loss parameter to achieve the same reflectance
# with a thinner PML, and if we use a too strong loss parameter, the discretized PML model
# deviates from the continuous model and we get stronger reflectance than the target value.
# Still, I found that the reflectance remains pretty low, and the convergence of iterative
# methods was more important.  So, I decided to specify the number of grid cells inside PML
# rather than the PML thickness.  This decision may change in the future.

export PMLParam
export get_pml_loc, gen_stretch_factor

mutable struct PMLParam
    m::Float  # degree of polynomial grading for σ and κ
    R::Float  # R: target reflection coefficient of PML
    κmax::Float  # maximum of real stretch factor
    amax::Float  # maximum of complex frequency shift
    ma::Float  # degree of polynomial grading for a
    PMLParam() = new(4.0, exp(-16), 1.0, 0.0, 4.0)
end

# Get the PML interface locations.
# The primal grid plane locations with the ghost plane locations are passed, from which the
# boundaries of the computational domain are retrieved.
function get_pml_loc(lg_prim::AbsVecReal, Npml::Tuple2{Integer})  # convenience method for K = 1
    lpml, Lpml = get_pml_loc((lg_prim,), SVec1.(Npml))

    return (lpml[nN][1],lpml[nP][1]), (Lpml[nN][1],Lpml[nP][1])
end
get_pml_loc(lg_prim::NTuple{K,AbsVecReal}, Npml::Tuple2{AbsVecInteger}) where {K} = get_pml_loc(lg_prim, SVector{K}.(Npml))

function get_pml_loc(lg_prim::NTuple{K,AbsVecReal},  # lg_prim[k] = primal vertex locations in k-direction (including ghost location at (+) end)
                     Npml::Tuple2{SVector{K,<:Integer}}) where {K}  # Npml[NEG][k] = number of cells inside PML at (-) end in k-direction
    bounds = (SVector{K}(getindex.(lg_prim,1)), SVector{K}(getindex.(lg_prim,lastindex.(lg_prim))))  # Tuple2{SVector{K,<:Real}}
    lprim = (l->@view(l[1:end-1])).(lg_prim)

    return get_pml_loc(lprim, bounds, Npml)
end

# Get the PML interface locations.
# The primal grid plane locations without the ghost plane locations are passed, so the
# positive-end boundaries of the computational domain cannot be retrieved from them.
# Therefore, the boundary information is passed separately.
function get_pml_loc(lprim::AbsVecReal, bounds::Tuple2{Real}, Npml::Tuple2{Integer})  # convenience method for K = 1
    lpml, Lpml = get_pml_loc((lprim,), SVec1.(bounds), SVec1.(Npml))

    return (lpml[nN][1],lpml[nP][1]), (Lpml[nN][1],Lpml[nP][1])
end
get_pml_loc(lprim::NTuple{K,AbsVecReal}, bounds::Tuple2{AbsVecReal}, Npml::Tuple2{AbsVecInteger}) where {K} =
    get_pml_loc(lprim, SVector{K}.(bounds), SVector{K}.(Npml))

function get_pml_loc(lprim::NTuple{K,AbsVecReal},  # lprim[k] = primal vertex locations in k-direction (excluding ghost location at (+) end)
                     bounds::Tuple2{SVector{K,<:Real}},  # bounds[NEG][k] = boundary of domain at (-) end in k-direction
                     Npml::Tuple2{SVector{K,<:Integer}}) where {K}  # Npml[NEG][k] = number of cells inside PML at (-) end in k-direction
    # Below, if Npml[nP] = 0, lprim[end+1-Npml[nP]] = lprim[end+1] is out of bound.  Because
    # lprim[end+1] points to the ghost point, we return the positive-end boundary in this
    # case.
    lpml = (map((v,n)->v[1+n], lprim, Npml[nN]), map((v,b,n)->(n==0 ? b : v[end+1-n]), lprim, bounds[nP], Npml[nP]))  # Tuple2{SVector{K,Float}}
    Lpml = (lpml[nN]-bounds[nN], bounds[nP]-lpml[nP])  # Tuple2{SVector{K,Float}}

    return lpml, Lpml
end

function calc_stretch_factor(ω::Number,  # angular frequency
                             d::Real,  # depth into PML
                             Lpml::Real,  # thickness of PML
                             pml::PMLParam  # PML parameters
                            )
    σmax = -(pml.m+1) * log(pml.R) /2Lpml  # -(m+1) ln(R) / (2 η Lpml), where η = 1 in units of η₀
    σ = σmax * (d/Lpml)^pml.m
    κ = 1 + (pml.κmax-1) * (d/Lpml)^pml.m
    a = pml.amax * (1 - d/Lpml)^pml.ma

    s_factor = κ + σ / (a + im*ω)  # s = κ + σ/(a + i ω ε), where ε = 1 in units of ε₀
end

function gen_stretch_factor(ω::Number, l::Tuple2{AbsVecReal}, lpml::Tuple2{Real}, Lpml::Tuple2{Real}, pml::PMLParam=PMLParam())
    sfactor = gen_stretch_factor(ω, ((l[nPR],), (l[nDL],)), SVec1.(lpml), SVec1.(Lpml), pml)

    return (sfactor[nPR][1], sfactor[nDL][1])
end

gen_stretch_factor(ω::Number, l::Tuple2{NTuple{K,AbsVecReal}}, lpml::Tuple2{AbsVecReal},
                   Lpml::Tuple2{AbsVecReal}, pml::PMLParam=PMLParam()) where {K} =
    gen_stretch_factor(ω, l, SVector{K}.(lpml), SVector{K}.(Lpml), pml)

function gen_stretch_factor(ω::Number,  # angular frequency
                            l::Tuple2{NTuple{K,AbsVecReal}},  # locations of primal and dual grid points
                            lpml::Tuple2{SVector{K,<:Real}},  # locations of PML interfaces
                            Lpml::Tuple2{SVector{K,<:Real}},  # thicknesses of PML
                            pml::PMLParam=PMLParam()  # PML parameters
                           ) where {K}
    N = length.(l[nPR])  # (Nx, Ny, Nz)
    sfactor = (ones.(CFloat,N), ones.(CFloat,N))
    for k = 1:K
        lpmlₙ, lpmlₚ = lpml[nN][k], lpml[nP][k]
        Lpmlₙ, Lpmlₚ = Lpml[nN][k], Lpml[nP][k]
        for ngt = nPD
            r = l[ngt][k]
            s = sfactor[ngt][k]
            for i = 1:N[k]
                r[i] < lpmlₙ && (s[i] = calc_stretch_factor(ω, lpmlₙ-r[i], Lpmlₙ, pml))
                r[i] > lpmlₚ && (s[i] = calc_stretch_factor(ω, r[i]-lpmlₚ, Lpmlₚ, pml))
            end
        end
    end

    return sfactor
end
