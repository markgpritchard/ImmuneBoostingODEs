
using DrWatson

@quickactivate :ImmuneBoostingODEs

using ModelingToolkit
using StructuralIdentifiability

@parameters β0, β1, ϕ, γ, μ, ψ, ω, proportiondetected
@independent_variables t
@variables S(t), I(t), R1(t), R2(t), x1(t), x2(t), y(t)
D = Differential(t)
β = β0 * (1 + β1 * x1)
λ = β * I
eqs = [
    D(S) ~ 3 * ω * (1 - S - I - R1 - R2) - λ * S + μ * (1 - S),
    D(I) ~ λ * S - (γ + μ) * I,
    D(R1) ~ γ * I + λ * ψ * (1 - S - I - R1) - (3 * ω + μ) * R1, 
    D(R2) ~ 3 * ω * R1 - (3 * ω + λ * ψ + μ) * R2,
    #D(R3) ~ 3 * ω * R2 - (3 * ω + λ * ψ + μ) * R3, 
    D(x1) ~ -2π * x2,
    D(x2) ~ 2π * x1,
]
measured_quantities = [y ~ λ * S * proportiondetected]
model = ODESystem(eqs, t; name=:SIRRRS)

assess_identifiability(model; measured_quantities)
