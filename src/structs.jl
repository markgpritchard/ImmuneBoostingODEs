
abstract type AbstractParameters end

@auto_hash_equals struct SirnsParameters{S, T, U, V} <: AbstractParameters where {S, T, U, V}
    β0::S
    β1::T
    ϕ::T
    γ::T
    μ::Float64 
    ψ::T
    ω::U
    originalβ0::V
    betaprimemultiplier::T
    finalbetaprime::T
    proportiondetected::T

    function SirnsParameters(
        β0::S, 
        β1::T, 
        ϕ::T, 
        γ::T, 
        μ, 
        ψ::T, 
        ω::U, 
        originalβ0::V, 
        betaprimemultiplier::T,
        finalbetaprime::T, 
        proportiondetected::T
    ) where {S, T, U, V}
        β0 >= 0 || throw(DomainError(β0, "β0 must not be negative"))
        0 <= β1 <= 1 || throw(_proportionerror(β1, "β1"))
        -π <= ϕ <= π || throw(DomainError(ϕ, "ϕ must be between -π and π"))
        γ >= 0 || throw(DomainError(γ, "γ must not be negative"))
        μ >= 0 || throw(DomainError(μ, "μ must not be negative"))
        ψ >= 0 || throw(DomainError(ψ, "ψ must not be negative"))
        ω >= 0 || throw(DomainError(ω, "ω must not be negative"))
        originalβ0  >= 0 || throw(DomainError(originalβ0, "originalβ0 must not be negative"))
        0 <= betaprimemultiplier <= 1 || throw(
            _proportionerror(betaprimemultiplier, "betaprimemultiplier")
        )
        finalbetaprime >= 0 || throw(
            DomainError(finalbetaprime, "finalbetaprime must not be negative")
        )
        0 <= proportiondetected <= 1 || throw(
            _proportionerror(proportiondetected, "proportiondetected")
        )
        return new{S, T, U, V}(
            β0, 
            β1, 
            ϕ, 
            γ, 
            μ, 
            ψ, 
            ω, 
            originalβ0, 
            betaprimemultiplier, 
            finalbetaprime, 
            proportiondetected
        )
    end
end     

function SirnsParameters(β0, γ::T, μ, ψ::T, ω) where T  
    β1 = zero(T)
    ϕ = zero(T)
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω) 
end 

function SirnsParameters(β0, β1::T, ϕ::T, γ::T, μ, ψ::T, ω) where T  
    betaprimemultiplier = one(T)
    finalbetaprime = one(T)
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, β0, betaprimemultiplier, finalbetaprime) 
end

function SirnsParameters(
    β0, 
    β1::T, 
    ϕ::T, 
    γ::T, 
    μ,
    ψ::T,
    ω, 
    originalβ0, 
    betaprimemultiplier::T, 
    finalbetaprime::T
) where T  
    proportiondetected = one(T)
    return SirnsParameters(
        β0, 
        β1, 
        ϕ, 
        γ, 
        μ, 
        ψ, 
        ω, 
        β0, 
        betaprimemultiplier, 
        finalbetaprime, 
        proportiondetected
    ) 
end

function SirnsParameters(
    β0, 
    β1::T, 
    ϕ::T, 
    γ::T, 
    μ, 
    ψ::Integer, 
    ω, 
    originalβ0, 
    betaprimemultiplier::T,
    finalbetaprime::T, 
    proportiondetected::T
) where T
    return SirnsParameters(
        β0, 
        β1, 
        ϕ, 
        γ, 
        μ, 
        T(ψ), 
        ω, 
        originalβ0, 
        betaprimemultiplier, 
        finalbetaprime, 
        proportiondetected
    )
end

function _proportionerror(parameter, parameterstring)
    _text = "$parameterstring is a proportion and must be between 0 and 1"
    return DomainError(parameter, _text)
end
