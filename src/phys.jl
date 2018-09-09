# Consider using Unitful.jl.

export PhysUnit, Oscillation  # types
export c₀, μ₀, ε₀, η₀, ℎ, ℏ  # constants
export in_L₀, in_ω₀, in_eV  # functions

const c₀ = 2.99792458e8  # in m/s
const μ₀ = 4e-7π  # in H/m
const ε₀ = 1 / (c₀^2 * μ₀)  # in F/m
const η₀ = √(μ₀ /ε₀)  # in Ohm
const ℎ = 4.13566733e-15  # in eV·sec
const ℏ = ℎ / 2π  # in eV·sec

struct PhysUnit
    L::Float
    ω::Float
    ε::Float
    μ::Float
    E::Float
    H::Float
    J::Float
    M::Float
    I::Float
    V::Float
    S::Float
    P::Float
    u::Float

    function PhysUnit(L₀::Real)
        L₀ > 0 || throw(ArgumentError("L₀ = $L₀ must be positive."));
        L = L₀;
        ω = c₀ / L;  # frequency in rad/s (L0-dependent)
        ε = ε₀;  # permittivity in eps0
        μ = μ₀;  # permeability in mu0
        E = 1.0;  # E-field in V/m
        H = E / η₀;  # H-field in A/m
        J = H / L;  # electric current density in A/m^2 (L0-dependent)
        M = E / L;  # magnetic current density in A/m^2 (L0-dependent)
        I = J * L^2;  # electric current in Amperes (L0-dependent)
        V = E * L;  # voltage in Volts (L0-dependent)
        S = E * H;  # Poynting vector in Watt/m^2
        P = S * L^2;  # power in Watt (L0-dependent)
        u = S / L;  # power density in Watt/m^3 (L0-dependent)

        new(L,ω,ε,μ,E,H,J,M,I,V,S,P,u)
    end
end

struct Oscillation
    λ::CFloat
    unit::PhysUnit
end

in_L₀(osc::Oscillation) = osc.λ
in_ω₀(osc::Oscillation) = 2π / osc.λ
in_eV(osc::Oscillation) = ℎ * c₀ / (osc.λ * osc.unit.L)
