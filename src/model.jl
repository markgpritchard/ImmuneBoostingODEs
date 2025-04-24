
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The ODE model 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


function sirns!(du, u, p, t)
    # Hard-coded to run with 3 resistant subcompartments 
    S, I, R1, R2, R3, x1, x2, = u 

    # transmission parameter
    β = p.β0 * (1 + p.β1 * x1) 
    λ = β * I
    
    du[1] = 3 * p.ω * R3 - λ * S + p.μ * (1 - S)                    # S
    du[2] = λ * S - (p.γ + p.μ) * I                                 # I
    du[3] = p.γ * I + λ * p.ψ * (R2 + R3) - (3 * p.ω + p.μ) * R1    # R1
    du[4] = 3 * p.ω * R1 - (3 * p.ω + λ * p.ψ + p.μ) * R2           # R2
    du[5] = 3 * p.ω * R2 - (3 * p.ω + λ * p.ψ + p.μ) * R3           # R3
    du[6] = -2π * x2                                                # x1
    du[7] = 2π * x1                                                 # x2
    du[8] = λ * S                                                   # cumulative cases 
end 

function _sirns!(du, u, p, t, λ)
    S, I, R1, R2, R3, x1, x2, cc = u
    
    du[1] = 3 * p.ω * R3 - λ * S + p.μ * (1 - S)                    # S
    du[2] = λ * S - (p.γ + p.μ) * I                                 # I
    du[3] = p.γ * I + λ * p.ψ * (R2 + R3) - (3 * p.ω + p.μ) * R1    # R1
    du[4] = 3 * p.ω * R1 - (3 * p.ω + λ * p.ψ + p.μ) * R2           # R2
    du[5] = 3 * p.ω * R2 - (3 * p.ω + λ * p.ψ + p.μ) * R3           # R3
    du[6] = -2π * x2                                                # x1
    du[7] = 2π * x1                                                 # x2
    du[8] = λ * S                                                   # cumulative cases 
end

sirns!(du, u, p::LambdaParms, t) = constantlambda_sirns!(du, u, p, t)
constantlambda_sirns!(du, u, p, t) = _sirns!(du, u, p, t, p.λ)


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

function sirns_u0(S0::S, I0, R1, R2, R3; p, t0 = 0) where S
    @assert +(S0, I0, R1, R2, R3) ≈ 1 "+($S0, $I0, $R1, $R2, $R3) = $(+(S0, I0, R1, R2, R3)) != 1"
    @assert min(S0, I0, R1, R2, R3) >= -1e-6 "min($S0, $I0, $R1, $R2, $R3) = $(min(S0, I0, R1, R2, R3)) < 0"
    u0 = Vector{S}(undef, 8)
    for (i, v) ∈ enumerate([ S0, I0, R1, R2, R3 ]) u0[i] = v end  
    u0[6] = cos(2π * t0 - p.ϕ)  # x1 
    u0[7] = sin(2π * t0 - p.ϕ)  # x2
    u0[8] = zero(S)  # cumulative cases
    return u0
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

function casespertimeblock(cc::Vector{T}) where T
    cases = Vector{T}(undef, length(cc) - 1)
    for i ∈ eachindex(cc)
        i == 1 && continue 
        newcases = cc[i] - cc[i-1] 
        # newcases should always be positive but occasionally the solver returns 
        # a very slightly negative value, such as -3e-298. Such values cannot be 
        # used with a Poisson distribution 
        if newcases < 0  
            cases[i-1] = zero(newcases) 
        else             
            cases[i-1] = newcases
        end
    end
    return cases
end


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

function reducetransmission!(integrator) 
    @unpack β0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0 = integrator.p
    integrator.p = SirnsParameters(
        reducedβ0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0
    )
end

# Callback function to restore βmean 

function restoretransmission!(integrator) 
    @unpack β0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0 = integrator.p
    integrator.p = SirnsParameters(
        restoredβ0, β1, ϕ, γ, μ, ψ, ω, originalβ0, reducedβ0, restoredβ0
    )
end
