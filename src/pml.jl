export PMLParam
export gen_stretch_factor

mutable struct PMLParam
    m::Float  # degree of polynomial grading for σ and κ
    R::Float  # R: target reflection coefficient of PML
    κmax::Float  # maximum of real stretch factor
    amax::Float  # maximum of complex frequency shift
    ma::Float  # degree of polynomial grading for a
    PMLParam() = new(4.0, exp(-16), 1.0, 0.0, 4.0)
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
