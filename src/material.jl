export Material
export matparam, kottke_avg_param
import Base:string

# Note that the output of tensorize is always a matrix: for K = 1, it is a 1×1 matrix.
tensorize(x::Number, ::Val{K}) where {K} = tensorize(@SVector(fill(complex(x), K)), Val(K))  # SComplex{3}(x,x,x): vector
tensorize(v::AbsVecNumber, ::Val{K}) where {K} = diagm(Val(0)=>SComplex{K}(v))  # complex(v) is allocating, so use SComplex{K}(v)
tensorize(m::AbsMatNumber, ::Val{K}) where {K} = SSComplex{K,K*K}(m)  # complex(m) is allocating, so use SSComplex{K,K*K}(m)

struct Material{Ke,Km,Ke²,Km²}  # e: electric, m: magnetic
    name::String
    param::Tuple{SSComplex{Ke,Ke²},SSComplex{Km,Km²}}  # (material parameter interacting with E, material parameter interacting with H)
    Material{Ke,Km,Ke²,Km²}(name, param) where {Ke,Km,Ke²,Km²} =
        new(name, (tensorize(param[1],Val(Ke)), tensorize(param[2],Val(Km))))  # suppress default outer constructor
end

Material{Ke,Km}(name::String, param::Tuple2{MatParam}) where {Ke,Km} = Material{Ke,Km,Ke*Ke,Km*Km}(name, param)
Material{Ke,Km}(name::String; ε::MatParam=1, μ::MatParam=1) where {Ke,Km} = Material{Ke,Km}(name, (ε,μ))  # ε, μ: keyword arguments
# Consider implementing Material2 and Material3 for 2D and 3D solvers.

string(m::Material) = m.name
matparam(m::Material, ft::FieldType) = m.param[Int(ft)]

kottke_avg_param(param1::AbsMatNumber, param2::AbsMatNumber, n12::AbsVecReal, rvol1::Real) =
    (K = length(n12); kottke_avg_param(SSComplex{K,K^2}(param1), SSComplex{K,K^2}(param2), SFloat{K}(n12), rvol1))

# Implement the averaging scheme between two local material parameter tensors developed in
# the paper by Kottke, Farjadpour, Johnson entitled "Perturbation theory for anisotropic
# dielectric interfaces, and application to subpixel smoothing of discretized numerical
# methods", Physical Review E 77 (2008): 036611.
function kottke_avg_param(param1::SSComplex{K}, param2::SSComplex{K}, n12::SFloat{K}, rvol1::Real) where {K}
    Scomp = @SMatrix rand(K,K-1)  # directions complementary to n12; works even for K = 1
    Stemp = [n12 Scomp]  # SMatrix{K,K}
    S = qr(Stemp).Q  # nonallocating; 1st column is normalized n12

    τ1 = τ_trans(transpose(S) * param1 * S)  # express param1 in S coordinates, and apply τ transform
    τ2 = τ_trans(transpose(S) * param2 * S)  # express param2 in S coordinates, and apply τ transform

    τavg = τ1 .* rvol1 + τ2 .* (1-rvol1)  # volume-weighted average

    return S * τ⁻¹_trans(τavg) * transpose(S)  # apply τ⁻¹ and transform back to global coordinates
end

# Kottke's averaging scheme reduced for the field dimension is orthogonal to the shape
# dimension (such as Hz in the TE mode).  In that case, the field is always parallel to the
# shape boundaries, so the averaging scheme reduces to simple arithmetic averaging.  For
# example for the case with Shape{2} and Kp = 1 (like Hz in the TE mode), the material
# parameter is an 1×1 tensor (= scalar) and therefore the situation is the same as averaging
# isotropic ε between E-field voxels, which is arithmetic averaging.  The only other cases
# are the case with Shape{1} and Kp = 2 (like Ex and Ey in the slaps in the z-direction) and
# the case with Shape{1} and Kp = 1  (like Ex in the slapes in the z-direction).  We can
# easily prove that in these cases again the averaging reduces to simple arithmetic
# averaging between material parameter tensors.  See Agenda > MaxwellFDM > Feb/22/2020.
kottke_avg_param(param1::AbsMatNumber, param2::AbsMatNumber, rvol1::Real) =
    param1 .* rvol1 + param2 .* (1-rvol1)

# Equation (4) of the paper by Kottke et al.  Need to verify for K = 1 and 2.
function τ_trans(ε::SSComplex3)
    ε₁₁, ε₂₁, ε₃₁, ε₁₂, ε₂₂, ε₃₂, ε₁₃, ε₂₃, ε₃₃ = ε

    return SSComplex3(
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

function τ_trans(ε::SSComplex2)
    ε₁₁, ε₂₁, ε₁₂, ε₂₂ = ε

    return SSComplex2(
        -1/ε₁₁, ε₂₁/ε₁₁,
        ε₁₂/ε₁₁, ε₂₂ - ε₂₁*ε₁₂/ε₁₁
    )
end

τ_trans(ε::SSComplex1) = -1 ./ ε

# function τ_trans(ε::SSComplex{K}) where {K}
#     τ_temp1 = @SMatrix [ε[i,j] for i = 2:K, j = 2:K]  # cannot use type parameter K here: error is generated
#     τ_temp2 = [@SMatrix(zeros(CFloat,1,K-1)); τ_temp1]
#     τ_temp = [@SVector(zeros(CFloat,K)) τ_temp2]
#
#     σcol_temp = @SVector [ε[i,1] for i = 2:K]
#     σcol = [SVector{1,CFloat}(-1); σcol_temp]
#
#     σrow_temp = @SVector [ε[1,j] for j = 2:K]
#     σrow = [SVector{1,CFloat}(-1); σrow_temp]
#
#     τ = τ_temp - (σcol * transpose(σrow)) ./ ε[1,1]
#
#     return τ
# end

# Equation (23) of the paper by Kottke et al.
function τ⁻¹_trans(τ::SSComplex3)
    τ₁₁, τ₂₁, τ₃₁, τ₁₂, τ₂₂, τ₃₂, τ₁₃, τ₂₃, τ₃₃ = τ

    return SSComplex3(
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

function τ⁻¹_trans(τ::SSComplex2)
    τ₁₁, τ₂₁, τ₁₂, τ₂₂ = τ

    return SSComplex3(
        -1/τ₁₁, -τ₂₁/τ₁₁
        -τ₁₂/τ₁₁, τ₂₂ - τ₂₁*τ₁₂/τ₁₁
    )
end

τ⁻¹_trans(τ::SSComplex1) = -1 ./ τ
