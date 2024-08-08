
abstract type AbstractParameters end

struct SirnsParameters <: AbstractParameters 
    β0          :: Float64
    β1          :: Float64
    ϕ           :: Float64
    γ           :: Float64
    μ           :: Float64 
    ψ           :: Float64
    ω           :: Float64
    βreduction1 :: Float64
    βreduction2 :: Float64
end     

struct LambdaParms <: AbstractParameters 
    λ           :: Float64 
    γ           :: Float64 
    μ           :: Float64 
    ψ           :: Float64 
    ω           :: Float64 
end   

mutable struct SirrrsParameters{S} <: AbstractParameters where S
    τ   :: S
    c   :: Matrix{Float64}
    γ   :: Float64
    ν   :: Float64
    μ   :: Vector{Float64} 
    ψ   :: S
    ω   :: Float64
end    


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Additional functions for structs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# When parameters are needed for a model with constant transmission parameter then it is 
# more convenient not to need to enter values for A and δ. Equally, if no non-pharmaceutical
# interventions are simulated, do not need to provide βreduction1 or βreduction2 explicitly

function SirnsParameters(β0, γ, μ, ψ, ω)  
    β1 = 0.0
    ϕ = 0.0
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω) 
end 

function SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω) 
    βreduction1 = 1.0
    βreduction2 = 1.0
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction1, βreduction2) 
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Imported functions for structs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function ==(a::SirrrsParameters, b::SirrrsParameters)
    return a.τ == b.τ && a.c == b.c && a.γ == b.γ && a.ν == b.ν && a.μ == b.μ && a.ψ == b.ψ && a.ω == b.ω 
end

function hash(a::SirrrsParameters)
    x = hash(a.τ)
    for i ∈ [ :c, :γ, :ν, :μ, :ψ, :ω ]
        x = hash(getproperty(a, i), x)
    end
    return x
end
