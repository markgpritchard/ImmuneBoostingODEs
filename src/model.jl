
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The ODE model 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function sirns!(du, u, p::AbstractParameters, t)
    S, I, R1, R2, R3, x1, x2, = u
    λ = _sirnslambda(p, u)
    
    du[1] = 3 * p.ω * R3 - λ * S + p.μ * (1 - S)  # S
    du[2] = λ * S - (p.γ + p.μ) * I  # I
    du[3] = p.γ * I + λ * p.ψ * (R2 + R3) - (3 * p.ω + p.μ) * R1  # R1
    du[4] = 3 * p.ω * R1 - (3 * p.ω + λ * p.ψ + p.μ) * R2  # R2
    du[5] = 3 * p.ω * R2 - (3 * p.ω + λ * p.ψ + p.μ) * R3  # R3
    du[6] = -2π * x2  # x1
    du[7] = 2π * x1  # x2
    du[8] = λ * S  # cumulative cases 
end 

function _sirnslambda(p, u)
    β = _sirnsbeta(p, u)
    I = u[2]
    λ = β * I
    return λ 
end

function _sirnsbeta(p::SirnsParameters, u)
    x1 = u[6]
    return _sirnsbeta(p.β0, p.β1, x1) 
end

function _sirnsbeta(p::AbstractVector, u)
    β0, β1, = p
    x1 = u[6]
    return _sirnsbeta(β0, β1, x1) 
end

_sirnsbeta(β0, β1, x1) = β0 * (1 + β1 * x1)

#@memoize function transformedsirns!(du, u, p, t)
function transformedsirns!(du, u, p, t)
    newparms = transformparameters(p)
    sirns!(du, u, newparms, t)
end

function transformparameters(p; r0, omega, gamma=48.7, mu=0.0087)
    r0 >= 0 || DomainError(r0, "r0 must not be negative")
    logitβ1, ϕ, logγ, logψ, logitβ′multiplier, logitfinalβ′, logitproportiondetected, β′, = p
    r0 * (exp(logγ) + mu) * β′ >= 0 || @warn "Negative β0 with p=$p, r0=$r0"
    return SirnsParameters(
        max(r0 * (exp(logγ) + mu) * β′, zero(r0 * (gamma + mu) * β′)),  # β0::T
        _logistic(logitβ1),  # β1::T
        ϕ,  # ϕ::T
        exp(logγ),  # γ::Float64
        mu,  # μ::Float64 
        exp(logψ),  # ψ::T
        omega,  # ω::T
        r0 * (exp(logγ) + mu),  # originalβ0::T
        _logistic(logitβ′multiplier),  # betaprimemultiplier::T
        _logistic(logitfinalβ′),  # finalbetaprime::T
        _logistic(logitproportiondetected),  # proportiondetected::T
    ) 
end

function transformparameterswithoutbetaprime(p; r0, omega, gamma=48.7, mu=0.0087)
    logitβ1, ϕ, logγ, logψ, logitβ′multiplier, logitfinalβ′, logitproportiondetected, = p
    return SirnsParameters(
        r0 * (exp(logγ) + mu),  # β0::T
        _logistic(logitβ1),  # β1::T
        ϕ,  # ϕ::T
        exp(logγ),  # γ::Float64
        mu,  # μ::Float64 
        exp(logψ),  # ψ::T
        omega,  # ω::T
        r0 * (exp(logγ) + mu),  # originalβ0::T
        _logistic(logitβ′multiplier),  # betaprimemultiplier::T
        _logistic(logitfinalβ′),  # finalbetaprime::T
        _logistic(logitproportiondetected),  # proportiondetected::T
    ) 
end

_logistic(x) = 1 / (1 + exp(-x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run models 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function run_sirns(
    u0::Vector{<:Real}, p::AbstractParameters, duration::Real; 
    t0=0, kwargs...
)
    tspan = ( Float64(t0), Float64(duration) )
    return run_sirns(u0, p, tspan; kwargs...)
end 

function run_sirns(
    u0::Vector{<:Real}, p::AbstractParameters, tspan::Tuple{<:Real, <:Real}; 
    kwargs...
)
    ts = ( Float64(tspan[1]), Float64(tspan[2]) )
    return run_sirns(u0, p, ts; kwargs...)
end 

function run_sirns(
    u0::Vector{<:Real}, p::AbstractParameters, tspan::Tuple{Float64, Float64}; 
    abstol=1e-12, alg=Vern9(lazy=false), maxiters=1e5, reltol=1e-12, saveat=0.0005, 
    kwargs...
)
    prob = ODEProblem(sirns!, u0, tspan, p)
    sol = solve(prob, alg; abstol, maxiters, reltol, saveat, kwargs...)
    return sol
end 

function sirns_u0(S0::S, I0::S; p, equalrs=false, kwargs...) where S
    Rtotal = 1 - (S0 + I0)
    if equalrs 
        rs = Rtotal / 3
        return sirns_u0(S0, I0, rs, rs, rs; p, kwargs...)
    else 
        return sirns_u0(S0, I0, Rtotal, zero(S), zero(S); p, kwargs...)
    end 
end 

function sirns_u0(S0, I0, R1, R2, R3; p, t0=0) 
    return _sirns_u0(S0, I0, R1, R2, R3, p, t0)
end 

function _sirns_u0(S0, I0, R1, R2, R3, p::SirnsParameters, t0)
    phi = p.ϕ
    return __sirns_u0(S0, I0, R1, R2, R3, phi, t0)
end 

function _sirns_u0(S0, I0, R1, R2, R3, p::SirnsParameters{T}, t0) where T <: ForwardDiff.Dual
    phi = ForwardDiff.value(p.ϕ)
    return __sirns_u0(S0, I0, R1, R2, R3, phi, t0)
end 

function __sirns_u0(S0::S, I0, R1, R2, R3, phi::Number, t0) where S
    @assert +(S0, I0, R1, R2, R3) ≈ 1 "+($S0, $I0, $R1, $R2, $R3) = $(+(S0, I0, R1, R2, R3)) != 1"
    @assert min(S0, I0, R1, R2, R3) >= -1e-6 "min($S0, $I0, $R1, $R2, $R3) = $(min(S0, I0, R1, R2, R3)) < 0"
    u0 = [
        S0,
        I0,
        R1,
        R2,
        R3,
        cos(2π * t0 - phi),  # x1 
        sin(2π * t0 - phi),  # x2
        zero(S),  # cumulative cases
    ]
    return u0
end 

function sirns_u0_transformedp(args...; p, r0, omega, kwargs...)
    newparms = transformparameterswithoutbetaprime(p; r0, omega)
    return sirns_u0(args...; p=newparms, kwargs...)
end

function modelcompartments(sol, p::T) where T <: Union{<:AbstractParameters, <:NamedTuple}
    inds = compartmentinds(sol)
    gt = sol.t[inds]
    S = modelcompartments(sol, 1, inds)
    I = modelcompartments(sol, 2, inds)
    R1 = modelcompartments(sol, 3, inds)
    R2 = modelcompartments(sol, 4, inds)
    R3 = modelcompartments(sol, 5, inds)
    Rtotal = @. R1 + R2 + R3
    β = p.β0 .* (1 .+ p.β1 .* modelcompartments(sol, 6, inds))
    cc = modelcompartments(sol, 8, inds)
    λ = β .* I
    return @dict gt S I R1 R2 R3 Rtotal cc β λ
end 

function modelcompartments(sol, i)
    inds = compartmentinds(sol)
    return modelcompartments(sol, i, inds)
end 

function modelcompartments(sol, c::Symbol, inds)
    i = COMPARTMENTINDICES[c]
    return modelcompartments(sol, i, inds)
end

modelcompartments(sol, i::Int, inds) = [ sol.u[j][i] for j ∈ inds ]
modelcompartments(sol, v::Vector{<:Integer}, inds) = [ sum(sol.u[j][v]) for j ∈ inds ]

function compartmentinds(sol)
    _gt = sol.t
    inds = findall(x -> x >= 0, _gt)
    return inds 
end 

casespertimeblock(d::Dict{Symbol, <:Any}) = casespertimeblock(d[:cc])
casespertimeblock(d::Dict{<:AbstractString, <:Any}) = casespertimeblock(d["cc"])
casespertimeblock(cc::Vector) = [ _newcases(cc, t) for t ∈ 2:length(cc) ]

_newcases(cc, t) = cc[t] - cc[t-1] < 0 ? zero(cc[t] - cc[t-1]) : cc[t] - cc[t-1]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to run the simulations with DrWatson.produce_or_load
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function modelincidence(
    p::AbstractParameters; 
    equalrs=true, I0=0.001, S0=0.5, tspan=( -1000.0, 10.0 )   
)
    u0 = sirns_u0(S0, I0; equalrs, p)
    sol = run_sirns(u0, p, tspan)
    cc = modelcompartments(sol, :cc)
    incidence = casespertimeblock(cc)
    gt = sol.t[compartmentinds(sol)]
    popfirst!(gt)
    return @ntuple gt incidence 
end

function pl_modelincidence(config::Dict{Symbol, <:Any})
    @unpack β0, β1, ϕ, γ, μ, ψ, ω, kw = config
    result = modelincidence(SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω); kw...)
    return tostringdict(result)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Callback functions 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Default callback functions used to simulate the effect of non-pharmaceutical interventions

reducetransmission!(integrator) = integrator.p = _reducetransmissionp(integrator.p)

function _reducetransmissionp(p::SirnsParameters)
    @unpack β0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0 = p
    newp = SirnsParameters(
        reducedβ0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0
    )
    return newp
end

function _reducetransmissionp(p::Tuple) 
    logr0, logitβ1, logitϕ, logψ, logitreduction1, logitreduction2, logitdetection = p
    newlogr0 = log(_logistic(logitreduction1)) + logr0
    return ( 
        newlogr0, 
        logitβ1, 
        logitϕ, 
        logψ, 
        logitreduction1, 
        logitreduction2, 
        logitdetection 
    )
end

function _reducetransmissionp(p::AbstractVector) 
    oldtuple = Tuple(p)
    newtuple = _reducetransmissionp(oldtuple)
    return [ newtuple... ]
end

# Callback function to restore βmean 

restoretransmission!(integrator) = integrator.p = _restoretransmissionp(integrator.p)

function _restoretransmissionp(p::SirnsParameters)
    @unpack β0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0 = p
    newp = SirnsParameters(
        restoredβ0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0
    )
    return newp
end

function _restoretransmissionp(p::Tuple) 
    logr0, logitβ1, logitϕ, logψ, logitreduction1, logitreduction2, logitdetection = p
    originallogr0 = logr0 - log(_logistic(logitreduction1))
    newlogr0 = log(_logistic(logitreduction2)) + originallogr0
    return ( 
        newlogr0, 
        logitβ1, 
        logitϕ, 
        logψ, 
        logitreduction1, 
        logitreduction2, 
        logitdetection 
    )
end

function _restoretransmissionp(p::AbstractVector) 
    oldtuple = Tuple(p)
    newtuple = _restoretransmissionp(oldtuple)
    return [ newtuple... ]
end
