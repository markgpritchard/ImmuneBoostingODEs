
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abstract type AbstractParameters end

struct SirnsParameters{T, U, V, W, X} <: AbstractParameters 
    β0::T
    β1::U
    ϕ::U
    γ::Float64
    μ::Float64 
    ψ::V
    ω::W
    βreductionfactor::T
    βreduction::X

    function SirnsParameters(
        β0::T, β1::U, ϕ, γ, μ, ψ::V, ω::W, βreductionfactor, βreduction
    ) where {T, U, V, W}
        for p in [β0, β1, γ, μ, ψ, ω, βreductionfactor, βreduction]  
            isnothing(p) && continue
            p >= 0 || throw(ArgumentError("$p, parameters must not be negative"))
        end
        # note, ϕ can be negative
        β1 <= 1 || throw(ArgumentError("$β1, β1 must not be greater than 1"))

        # type of X has to be able to include values reached after the callback  
        X = typeof(exp(log(βreduction) * βreductionfactor))

        return new{T, U, V, W, X}(
            β0, 
            β1, 
            convert(U, ϕ), 
            convert(Float64, γ), 
            convert(Float64, μ), 
            ψ, 
            ω, 
            convert(T, βreductionfactor),
            convert(X, βreduction),
        )
    end

    function SirnsParameters(
        p::SirnsParameters{T, U, V, W, X}, βreduction
    ) where {T, U, V, W, X}
        # internal function without checks on parameters that have already been checked
        newbetareduction = exp(log(βreduction) * p.βreductionfactor)

        return new{T, U, V, W, typeof(newbetareduction)}(p.β0, p.β1, p.ϕ, p.γ, p.μ, p.ψ, p.ω, p.βreductionfactor, βreduction)
    end
end 

function SirnsParameters( ; 
    β0=0, β1=nothing, ϕ=nothing, γ=0, μ=0, ψ=0, ω=0, βreductionfactor=1, βreduction=1,
)
    return SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreductionfactor, βreduction)
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
    a.βreduction == b.βreduction || return false 
    return true 
end

function Base.hash(a::SirnsParameters, h)
    h = hash(:SirnsParameters, h)
    for i ∈ [:β0, :β1, :ϕ, :γ, :μ, :ψ, :ω, :βreduction]
        h = hash(getproperty(a, i), h)
    end
    return h
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The ODE model 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    λ = _sirnslambda(p, u)
    return _sirns!(du, u, p, t, λ)
end 

function _sirns!(du, u, p, t, λ)
    # Hard-coded to run with 3 resistant subcompartments 
    S, I, R1, R2, R3, x1, x2, cc = u
    
    du[1] = 3 * p.ω * R3 - λ * S + p.μ * (1 - S)  # S
    du[2] = λ * S - (p.γ + p.μ) * I  # I
    du[3] = p.γ * I + λ * p.ψ * (R2 + R3) - (3 * p.ω + p.μ) * R1  # R1
    du[4] = 3 * p.ω * R1 - (3 * p.ω + λ * p.ψ + p.μ) * R2  # R2
    du[5] = 3 * p.ω * R2 - (3 * p.ω + λ * p.ψ + p.μ) * R3  # R3
    du[6] = -2π * x2  # x1
    du[7] = 2π * x1  # x2
    du[8] = λ * S  # cumulative cases
    return nothing 
end

sirns!(du, u, p::LambdaParms, t) = _sirns!(du, u, p, t, p.λ)

"""
    constantlambda_sirns!(du, u, p, t)

`p` is expected to be of type `LambdaParms` but other types will be accepted provided 
    that all required parameters are named.

The function `sirns_u0` can be used to produce an appropriate vector for `u`. Note 
    that this vector must include all model compartments, the transmission parameter, 
    and a value of cumulative infections.

See also `sirns!`.
"""
constantlambda_sirns!(du, u, p, t) = _sirns!(du, u, p, t, p.λ)
# Hard-coded to run with 3 resistant subcompartments 

_sirnsbeta_0(p) = p.β0 * p.βreduction
_sirnsbeta_0(p::SirnsParameters{<:Any, <:Any, <:Any, <:Any, Nothing}) = p.β0
_sirnsbeta(p, u) = _sirnsbeta_0(p) * (1 + p.β1 * u[6])  # β0 * (1 + β1 * x1)

function _sirnsbeta(p::SirnsParameters{<:Any, Nothing, <:Any, <:Any, <:Any}, ::Any)
    return _sirnsbeta_0(p)
end

_sirnslambda(p, u) = _sirnsbeta(p, u) * u[2]  # λ = β * I


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
function run_sirns(u0::Vector{<:Real}, p::AbstractParameters, duration::Real; t0=0, kwargs...)
    return run_sirns(u0, p, (t0, duration); kwargs...)
end 

function run_sirns(
    u0::Vector{<:Real}, p::AbstractParameters, tspan::Tuple; 
    abstol=1e-12, alg=Vern9(; lazy=false), maxiters=1e5, reltol=1e-12, saveat=0.0005, 
    kwargs...
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
function sirns_u0(S0::S, I0::S; p, equalrs=false, kwargs...) where S
    Rtotal = 1 - (S0 + I0)
    if equalrs 
        rs = Rtotal / 3
        return sirns_u0(S0, I0, rs, rs, rs; p, kwargs...)
    else 
        return sirns_u0(S0, I0, Rtotal, zero(S), zero(S); p, kwargs...)
    end 
end 

function sirns_u0(S0::S, I0, R1, R2, R3; p, t0=0) where S
    #= if isnan(S0)
        @warn "NaN values passed to u0, S0=$S0, I0=$I0, R1=$R1, R2=$R2, R3=$R3, p=$p" 
    end =#
  #=  @assert +(S0, I0, R1, R2, R3) ≈ 1 "+($S0, $I0, $R1, $R2, $R3) = $(+(S0, I0, R1, R2, R3)) != 1"
    @assert min(S0, I0, R1, R2, R3) >= -1e-6 "min($S0, $I0, $R1, $R2, $R3) = $(min(S0, I0, R1, R2, R3)) < 0" =#
    +(S0, I0, R1, R2, R3) ≈ 1 || @warn "model expects compartment values to sum to 1: S0=$S0, I0=$I0, R1=$R1, R2=$R2, R3=$R3"
    min(S0, I0, R1, R2, R3) >= -1e-6 || @warn "model expects all compartment values to be ≥ 0: S0=$S0, I0=$I0, R1=$R1, R2=$R2, R3=$R3"
    u0 = [
        S0, 
        I0, 
        R1, 
        R2, 
        R3,
        _initialx1(S, p, t0),
        _initialx2(S, p, t0),
        zero(S)  # cumulative cases
    ]
    return u0
end 

_initialx1(::Any, p, t0) = cos(2π * t0 - p.ϕ)
_initialx1(S, ::SirnsParameters{<:Any, Nothing, <:Any, <:Any}, ::Any) = one(S) 
_initialx2(::Any, p, t0) = sin(2π * t0 - p.ϕ)
_initialx2(S, ::SirnsParameters{<:Any, Nothing, <:Any, <:Any}, ::Any) = zero(S) 

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

modelcompartments(sol, i::Int, inds) = [sol.u[j][i] for j ∈ inds]
modelcompartments(sol, v::Vector{<:Integer}, inds) = [sum(sol.u[j][v]) for j ∈ inds]

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
    equalrs=true, I0=0.001, S0=0.5, tspan=(-1000.0, 10.0)   
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
    result = modelincidence(SirnsParameters(; β0, β1, ϕ, γ, μ, ψ, ω); kw...)
    return tostringdict(result)
end
