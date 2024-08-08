
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The ODE model 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# New way of coding sirns! saves approximately 90% of the runtime. Tests added to 
# check that x1 == cos(2π t - ϕ)

""" 
    sirns!(du, u, p, t)

The ordinary differential equations model.

`p` is expected to be of type `SirnsParameters` or `LambdaParms` to run the model 
    with a constant force of infection. Other types will be accepted provided that 
    all required parameters are named.

The function `sirns_u0` can be used to produce an appropriate vector for `u`. Note 
    that this vector must include all model compartments, the transmission parameter, 
    and a value of cumulative infections.
"""
function sirns!(du, u, p, t)
    # Hard-coded to run with 3 resistant subcompartments 
    S, I, R1, R2, R3, x1, x2, cc = u 
    #@assert minimum(u[1:5]) >= 0

    # transmission parameter
    β = p.β0 * (1 + p.β1 * x1) # more efficient version than cos(2π(t-ϕ))
    @assert β >= 0 "β<0 when p=$p, x1=$x1"
    λ = β * I
    _sirns!(du, u, p, t, λ)
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

"""
    constantlambda_sirns!(du, u, p, t)

`p` is expected to be of type `LambdaParms` but other types will be accepted provided 
    that all required parameters are named.

The function `sirns_u0` can be used to produce an appropriate vector for `u`. Note 
    that this vector must include all model compartments, the transmission parameter, 
    and a value of cumulative infections.

See also `sirns!`.
"""
function constantlambda_sirns!(du, u, p, t)
    # Hard-coded to run with 3 resistant subcompartments 
    S, I, R1, R2, R3, x1, x2, cc = u 

    _sirns!(du, u, p, t, p.λ)
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run models 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""" 
    run_sirns(u0::Vector{<:Real}, p::AbstractParameters, duration::Real; t0 = 0, <keyword arguments>)
    run_sirns(u0::Vector{<:Real}, p::AbstractParameters, tspan::Tuple{Float64, Float64}; <keyword arguments>)

Run the `sirns!` model.

The function `sirns_u0` can be used to produce an appropriate vector for `u0`. Note 
    that this vector must include all model compartments, the transmission parameter, 
    and a value of cumulative infections.

## Keyword arguments

All keyword arguments are passed to `DifferentialEquations.solve`. The following 
    defaults are provided to reduce the likelihood that differences between simulations 
    are due to differences in running the ODE solver.

* `abstol = 1e-12` 
* `alg = Vern9(lazy = false)` (note this is a positional argument for `DifferentialEquations.solve`)
* `maxiters = 1e5`
* `reltol = 1e-12` 
* `saveat = .0005`
""" 
function run_sirns(u0::Vector{<:Real}, p::AbstractParameters, duration::Real; t0 = 0, kwargs...)
    tspan = ( Float64(t0), Float64(duration) )
    return run_sirns(u0, p, tspan; kwargs...)
end 

function run_sirns(u0::Vector{<:Real}, p::AbstractParameters, tspan::Tuple{<:Real, <:Real}; kwargs...)
    ts = ( Float64(tspan[1]), Float64(tspan[2]) )
    return run_sirns(u0, p, ts; kwargs...)
end 

function run_sirns(u0::Vector{<:Real}, p::AbstractParameters, tspan::Tuple{Float64, Float64}; 
        abstol = 1e-12, alg = Vern9(lazy = false), maxiters = 1e5, reltol = 1e-12, 
        saveat = .0005, kwargs...
    )
    prob = ODEProblem(sirns!, u0, tspan, p)
    sol = solve(prob, alg; abstol, maxiters, reltol, saveat, kwargs...)
    return sol
end 

"""
    sirns_u0(S0, I0[, Rs0...]; p, <keyword arguments>)

Construct the vector `u0` to run the model `sirns!`.

Inputs can be values for all compartments or just for `S0` and `I0`. If only `S0` 
    and `I0` are supplied, the additional keyword argument `equalrs = false` determines 
    whether the resistant subcompartments will take equal numbers (the default is 
    to put all resistant individuals into the first subcompartment).

## Keyword arguments

* `t0 = 0` to set the start time for the model so that initial values of `β` are calculated 
    correctly.
"""
function sirns_u0(S0::S, I0::S; p, equalrs = false, kwargs...) where S
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
    u0[8] = zero(S)               # cumulative cases
    return u0
end 

""" 
    modelcompartments(sol, p)
    modelcompartments(sol, c)
    modelcompartments(sol, c, inds)

Return vectors of compartment sizes from the ODE solver outputs. 

## Accepted inputs
* `sol` is the output solution from the ODE solver
* `p` is an `AbstractParameters` struct or a `NamedTuple` of model parameters. If 
    `p` is supplied, the output is a `Dict` of all model compartments, the total 
    proportions immune, and vectors of the transmission parameter and force of infection.
* `c` is either a symbol, an integer, or a vector of symbols or integers, indicating 
    which compartment's values to return. 
* `inds` is a vector of saved times to be included in the vector. The default is 
    all `t ≥ 0`.
"""
function modelcompartments(sol, p::T) where T <: Union{<:AbstractParameters, <:NamedTuple}
    inds    = compartmentinds(sol)
    gt      = sol.t[inds]
    S       = modelcompartments(sol, 1, inds)
    I       = modelcompartments(sol, 2, inds)
    R1      = modelcompartments(sol, 3, inds)
    R2      = modelcompartments(sol, 4, inds)
    R3      = modelcompartments(sol, 5, inds)
    Rtotal  = @. R1 + R2 + R3
    β       = p.β0 .* (1 .+ p.β1 .* modelcompartments(sol, 6, inds))
    cc      = modelcompartments(sol, 8, inds)
    λ       = β .* I
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

"""
    casespertimeblock(cc)

Calculates incidence from a vector of cumulative infection data.

The returned vector will be 1 shorter than the supplied vector as no incidence is 
    calculated for the initial time point. 

Input can be a vector of values, or a Dict containing an element labelled `:cc` or 
    `\"cc\"`.
"""
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

"""
    pl_modelincidence(config::Dict{Symbol, <:Any})

Function to calculate model incidence using `DrWatson.produce_or_load`.

`config` is a `Dict{Symbol, <:Any}` containing `β0`, `β1`, `ϕ`, `γ`, `μ`, `ψ`, `ω`, 
    which are passed to `SirnsParameters`, and `kw`, which is a `NamedTuple` of 
    keyword arguments. The default keyword arguments are, 
* `equalrs = true`, passed to `sirns_u0`
* `I0 = .001`, passed to `sirns_u0`
* `S0 = .5`, passed to `sirns_u0`
* `tspan = ( -1000., 10. )`, for the duration of the simulation 
"""
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
