export Material
export matparam, kottke_avg_param
import Base:string

tensorize(x::Number) = tensorize(SVec3Complex(x,x,x))
tensorize(v::AbsVecNumber) = diagm(Val(0)=>SVec3Complex(v))
tensorize(m::AbsMatNumber) = SMat3Complex(m)

struct Material
    name::String
    param::Tuple2{SMat3Complex}  # (material parameter interacting with E, material parameter interacting with H)
end
Material(name::String; ε::MatParam=1, μ::MatParam=1) = Material(name, (tensorize(ε), tensorize(μ)))

string(m::Material) = m.name
matparam(m::Material, ft::FieldType) = m.param[Int(ft)]

kottke_avg_param(param1::AbsMatNumber, param2::AbsMatNumber, n12::AbsVecReal, rvol1::Real) =
    kottke_avg_param(SMat3Complex(param1), SMat3Complex(param2), SVec3Float(n12), rvol1)

# Implement the averaging scheme between two local material parameter tensors developed in
# the paper by Kottke, Farjadpour, Johnson entitled "Perturbation theory for anisotropic
# dielectric interfaces, and application to subpixel smoothing of discretized numerical
# methods", Physical Review E 77 (2008): 036611.
function kottke_avg_param(param1::SMat3Complex, param2::SMat3Complex, n12::SVec3Float, rvol1::Real)
    n = normalize(n12)

    # Pick a vector that is not along n.
    if any(n .== 0)
    	h = (n .== 0)
    else
    	h = SVector(1., 0., 0.)
    end

    # Create two vectors that are normal to n and normal to each other.
    h = normalize(n×h)
    v = normalize(n×h)

    # Create a local Cartesian coordinate system.
    S = [n h v]  # unitary

    τ1 = τ_trans(transpose(S) * param1 * S)  # express param1 in S coordinates, and apply τ transform
    τ2 = τ_trans(transpose(S) * param2 * S)  # express param2 in S coordinates, and apply τ transform

    τavg = τ1 .* rvol1 + τ2 .* (1-rvol1)  # volume-weighted average

    return S * τ⁻¹_trans(τavg) * transpose(S)  # apply τ⁻¹ and transform back to global coordinates
end

# Equation (4) of the paper by Kottke et al.
function τ_trans(ε)
    ε₁₁, ε₂₁, ε₃₁, ε₁₂, ε₂₂, ε₃₂, ε₁₃, ε₂₃, ε₃₃ = ε

    return SMat3Complex(
        -1/ε₁₁, ε₂₁/ε₁₁, ε₃₁/ε₁₁,
        ε₁₂/ε₁₁, ε₂₂ - ε₂₁*ε₁₂/ε₁₁, ε₃₂ - ε₃₁*ε₁₂/ε₁₁,
        ε₁₃/ε₁₁, ε₂₃ - ε₂₁*ε₁₃/ε₁₁, ε₃₃ - ε₃₁*ε₁₃/ε₁₁
    )
    # t = zeros(3,3)
    # t[2:3, 2:3] = s[2:3, 2:3]
    # s11 = s[1,1]
    # s[1,1] = -1
    #
    # t = t - (s(:,1) * s(1,:)) ./ s11
    # return t
end

# Equation (23) of the paper by Kottke et al.
function τ⁻¹_trans(τ)
    τ₁₁, τ₂₁, τ₃₁, τ₁₂, τ₂₂, τ₃₂, τ₁₃, τ₂₃, τ₃₃ = τ

    return SMat3Complex(
        -1/τ₁₁, -τ₂₁/τ₁₁, -τ₃₁/τ₁₁,
        -τ₁₂/τ₁₁, τ₂₂ - τ₂₁*τ₁₂/τ₁₁, τ₃₂ - τ₃₁*τ₁₂/τ₁₁,
        -τ₁₃/τ₁₁, τ₂₃ - τ₂₁*τ₁₃/τ₁₁, τ₃₃ - τ₃₁*τ₁₃/τ₁₁
    )

    # s = zeros(3,3)
    # s[2:3, 2:3] = t[2:3, 2:3]
    # t11 = t[1,1]
    # t[1,1] = 1
    #
    # s = s - (t(:,1) * t(1,:)) ./ t11
    # return s
end
