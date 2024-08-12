
abstract type AbstractParameters end

struct SirnsParameters{T} <: AbstractParameters where T
    β0          :: T
    β1          :: T
    ϕ           :: T
    γ           :: Float64
    μ           :: Float64 
    ψ           :: T
    ω           :: Float64
    originalβ0  :: T
    reducedβ0   :: T
    restoredβ0  :: T
end     

struct LambdaParms <: AbstractParameters 
    λ           :: Float64 
    γ           :: Float64 
    μ           :: Float64 
    ψ           :: Float64 
    ω           :: Float64 
end   


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Additional functions for structs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# When parameters are needed for a model with constant transmission parameter then it is
# more convenient not to need to enter values for A and δ. Equally, if no non-pharmaceutical
# interventions are simulated, do not need to provide reducedβ0 or restoredβ0 explicitly

function SirnsParameters(β0::T, γ, μ, ψ::T, ω::T) where T  
    β1 = zero(T)
    ϕ = zero(T)
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω) 
end 

function SirnsParameters(β0::T, γ::T, μ::T, ψ::Integer, ω::T) where T <: Float64
    return SirnsParameters(β0, γ, μ, Float64(ψ), ω)
end

function SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω) 
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, β0, β0, β0) 
end

# function to ensure consistency with older code 
function SirnsParameters(β0::T, β1::T, ϕ::T, γ, μ, ψ::T, ω::T, βreduction1, βreduction2) where T
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, β0, β0 * βreduction1, β0 * βreduction2) 
end

function SirnsParameters(
    β0::T, β1::T, ϕ::T, γ, μ, ψ::Integer, ω::T, βreduction1, βreduction2
) where T <: Float64
    return SirnsParameters(β0, β1, ϕ, γ, μ, Float64(ψ), ω, βreduction1, βreduction2) 
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Imported functions for structs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to show that two sets of parameters are equal. May assist in identifying repeat
# instances for `@memoize`.

function ==(a::SirnsParameters, b::SirnsParameters)
    return a.β0 == b.β0 && 
        a.β1 == b.β1 && 
        a.ϕ == b.ϕ && 
        a.γ == b.γ && 
        a.μ == b.μ && 
        a.ψ == b.ψ && 
        a.ω == b.ω && 
        a.originalβ0 == b.originalβ0 && 
        a.reducedβ0 == b.reducedβ0 && 
        a.restoredβ0 == b.restoredβ0 
end

function hash(a::SirnsParameters)
    x = hash(a.β0)
    for i ∈ [ :β1, :ϕ, :γ, :μ, :ψ, :ω, :originalβ0, :reducedβ0, :restoredβ0 ]
        x = hash(getproperty(a, i), x)
    end
    return x
end
