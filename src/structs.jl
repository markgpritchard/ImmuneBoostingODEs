
abstract type AbstractParameters end

struct SirnsParameters{T} <: AbstractParameters 
    β0::T
    β1::T
    ϕ::T
    γ::Float64
    μ::Float64 
    ψ::T
    ω::Float64
    originalβ0::T
    reducedβ0::T
    restoredβ0::T

    function SirnsParameters(
        β0::T, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0
    ) where T
        for p in [β0, β1, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0]  # ϕ can be negative
            p >= 0 || throw(ArgumentError("$p, parameters must not be negative"))
        end
        return new{T}(
            β0, 
            convert(T, β1), 
            convert(T, ϕ), 
            convert(Float64, γ), 
            convert(Float64, μ), 
            convert(T, ψ), 
            convert(Float64, ω), 
            convert(T, originalβ0), 
            convert(T, reducedβ0), 
            convert(T, restoredβ0)
        )
    end
end 

function SirnsParameters(; 
    β0, β1=0, ϕ=0, γ=0, μ=0, ψ=0, ω=0, originalβ0=β0, reducedβ0=β0, restoredβ0=β0,
)
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0)
end

struct LambdaParms <: AbstractParameters 
    λ::Float64 
    γ::Float64 
    μ::Float64 
    ψ::Float64 
    ω::Float64 

    function LambdaParms(λ, γ, μ, ψ, ω)
        for p in [λ, γ, μ, ψ, ω]
            p >= 0 || throw(ArgumentError("$p, parameters must not be negative"))
        end
        return new(
            convert(Float64, λ), 
            convert(Float64, γ), 
            convert(Float64, μ), 
            convert(Float64, ψ), 
            convert(Float64, ω),
        )
    end
end   

function Base.:(==)(a::SirnsParameters, b::SirnsParameters)
    a.β0 == b.β0 || return false 
    a.β1 == b.β1 || return false 
    a.ϕ == b.ϕ || return false 
    a.γ == b.γ || return false 
    a.μ == b.μ || return false 
    a.ψ == b.ψ || return false 
    a.ω == b.ω || return false 
    a.originalβ0 == b.originalβ0 || return false 
    a.reducedβ0 == b.reducedβ0 || return false 
    a.restoredβ0 == b.restoredβ0 || return false 
    return true 
end

function Base.hash(a::SirnsParameters, h)
    x = hash(:SirnsParameters, h)
    for i ∈ [:β0, :β1, :ϕ, :γ, :μ, :ψ, :ω, :originalβ0, :reducedβ0, :restoredβ0]
        x = hash(getproperty(a, i), x)
    end
    return x
end
