
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Values at equilibrium
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
    equil(p; warntol = 1e-10)

Calculate proportions in each compartment at equilibrium.

Will return endemic equilibrium if able and disease-free equilibrium otherwise.

`p` is a `SirnsParameters` structure or other named list of parameters.

Function checks that the compartments sum to 1. If the difference from 1 is greater 
    than `warntol` a warning is displayed.
"""
function equil(p; warntol = 1e-10)
    u = zeros(5)
    I = equili(p)
    _equil!(u, I, p)
    equilwarning(u, warntol)
    return u
end 

function _equil!(u, I, p)
    if 0 < I <= 1 _endemicequil!(u, I, p)   # endemic equilibrium
    else          u[1] = 1.                 # disease-free equilibium  
    end 
end 

function _endemicequil!(u, I, p)
    u[1] = equils(p) 
    u[2] = I 
    for i ∈ 1:3 u[2+i] = equilri(p, I, i) end 
end 

# Provide a warning if needed by `equil`
function equilwarning(u, warntol)
    if sum(u) < 1 - warntol || sum(u) > 1 + warntol
        @warn """
            equil($p) -> [$u]; sum(u) = $(sum(u))
            `u` is still returned
        """
    end
end 

"""
    equils(p)
    equils(βmean, γ, μ)

Calculate the endemic equilibrium proportion susceptible.

Can accept `p` as a `SirnsParameters` structure or other named list of parameters, 
    or individual parameters `βmean`, `γ` and `μ`.

If there is no endemic equilibrium, it will return a value of `1` for the disease-free 
    equilibrium.
"""
equils(p) = equils(p.β0, p.γ, p.μ)

function equils(β0, γ, μ) 
    S = _equils(β0, γ, μ) 
    if S > 1 return one(S) 
    else     return S 
    end 
end 

_equils(β0, γ, μ) = (γ + μ) / β0

""" 
    equili(p) 

Calculate the proportion infectious at endemic equilibrium.

`p` is a SirnsParameters or other named list of parameters.

If there is no endemic equilibrium, it will return a value of `0` for the disease-free 
    equilibrium.
"""
function equili(p) 
    s1 = equiliproblem(0, p)
    s2 = equiliproblem(1, p)
    if s1 * s2 < 0 
        Z = ZeroProblem(equiliproblem, ( 0, 1 ))
        return solve(Z, p)
    else # s1 and s2 are the same sign then I* ∉ [0, 1] so no endemic equilibrium
        return .0 
    end
end 

## equilrn (not exported)
# Calculates the proportion of the population that would be in the final resistant 
# subcompartment given parameters, `p`, and a proportion infection, `I`.
equilrn(p, I) = equilrn(p.β0, p.γ, p.μ, p.ω, I)
equilrn(β0, γ, μ, ω, I) = ((I + μ / β0) * (γ + μ) - μ) / (3 * ω)

## equilri_multiplier
# Calculate the ratio of the proportion in the `i`th resistant subcompartment compared 
# to the final subcompartment at endemic equilibrium, given parameters `p` and a 
# proportion `I` in the infectious compartment.
equilri_multiplier(p, I, i) = equilri_multiplier(p.β0, p.μ, p.ψ, p.ω, I, i)

function equilri_multiplier(β0, μ, ψ, ω, I, i)
    power = 3 - i 
    return ((β0 * ψ * I + 3 * ω + μ) / (3 * ω))^power 
end 

""" 
    equilri(p, I, i)

Calculate the proportion in the `i`th resistant subcompartment at endemic equilibrium, 
    given parameters `p` and a proportion `I` infectious.

This function may return values >0 even if there is no endemic equilibrium.
"""
equilri(p, I, i) = equilrn(p, I) * equilri_multiplier(p, I, i)

"""
    equilr(p, I)

Calculate the total proportion in the resistant subcompartments at endemic equilibrium, 
    given parameters `p` and a proportion `I` in the infectious compartment.

This function may return values >0 even if there is no endemic equilibrium.
"""
equilr(p, I) = sum([ equilri(p, I, i) for i ∈ 1:3 ])

## equiliproblem (not exported)
# When x is the equilibium proportion infectious for parameters `p` then this function 
# equals 0.
equiliproblem(x, p) = equils(p) + x + equilr(p, x) - 1


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Eigenvalues at equilibrium 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
    equileigen(p; sortby = λ -> (real(λ),imag(λ)), warntol = 1e-10)

Calculate eigenvalues of the Jacobian at equilibium.

`sortby` orders the output to aid plotting results. `warntol` is passed to `equil`.
"""
function equileigen(p; sortby = λ -> (real(λ),imag(λ)), warntol = 1e-10)
    u = equil(p; warntol)
    J = modeljacobian(p, u)
    return eigen(J; sortby).values 
end 

""" 
    modeljacobian(p, u)

Return the Jacobian matrix for the model with parameters `p` and compartment 
    values `u`
"""
function modeljacobian(p, u)
    J = zeros(5, 5)
    sirns_jac!(J, u, p)
    return J
end 

# This function calculates the Jacobian for a given set of parameters to allow calculation 
# of eigenvalues and understanding of stability of equilibria. Versions of the function 
# that had been designed to calculate the Jacobian while simulations were running 
# did not speed up the simulation and added a potential for introducing error to 
# the simulation so those have been deleted.
function sirns_jac!(J, u, p)
    sirns_jacwarning(p)
    S, I, = u
    β = p.β0
    # with repect to S
    J[1, 1] = -β * I - p.μ
    J[2, 1] = β * I
    # with respect to I 
    J[1, 2] = -β * S
    J[2, 2] = β * S - p.γ - p.μ
    J[3, 2] = p.γ + β * p.ψ * sum(u[4:5])
    for i ∈ 4:5 J[i, 2] = -β * p.ψ * u[i] end
    # with respect to R1 
    J[3, 3] = -3 * p.ω - p.μ
    J[4, 3] = 3 * p.ω 
    # with respect to further resistant subcompartments 
    for i ∈ 4:5
        J[3, i] = β * p.ψ * I 
        J[i, i] = -3 * p.ω - β * p.ψ * I - p.μ 
    end 
    J[5, 4] = 3 * p.ω  
    J[1, 5] = 3 * p.ω 
end 

function sirns_jacwarning(p)
    if p.β1 != 0
        @warn "sirns_jac! received p.β1 = $(p.β1) but will assume no seasonal forcing"
    end
end

""" 
    maxequileigen(p; <keyword arguments>)

Return the eigenvalue from the model's Jacobian with the largest real part.

All keyword arguments are passed to `equileigen`.
"""
function maxequileigen(p; kwargs...)
    eigens = equileigen(p; kwargs...)
    re = real.(eigens)
    ind = findmax(re)[2] # findmax returns a tuple of the maximum value and its index
    return eigens[ind]
end 

""" 
    realmaxequileigen(β, γ, μ, ψ, ω; n = 3, <keyword arguments>)
    realmaxequileigen(p; <keyword arguments>)

Return the real part of the eigenvalue from the model's Jacobian with the largest 
    real part.

All keyword arguments except `n` are passed to `maxequileigen`.
"""
realmaxequileigen(p; kwargs...) = real(maxequileigen(p; kwargs...))

function realmaxequileigen(β, γ, μ, ψ, ω; kwargs...)
    p = SirnsParameters(β, γ, μ, ψ, ω)
    return realmaxequileigen(p; kwargs...)
end
#=
""" 
    equilplotdata(p::SirnsParameters; <keyword arguments>)

Return data to plot whether endemic equilibrium is stable.

Returns:
* `-1` if no endemic equilibrium 
* `0` if the endemic equilibrium is unstable 
* `1` if the endemic equilibrium is stable

All keyword arguments are passed to `realmaxequileigen`.
"""
function equilplotdata(p::SirnsParameters; kwargs...)
    I = equili(p) 
    if I == 0 
        return -1               # no endemic equilibium
    else 
        eig = realmaxequileigen(p; kwargs...)
        if eig < 0  return 1    # stable endemic equilibium
        else        return 0    # unstable endemic equilibium
        end
    end
end
=#
""" 
    findpsi(β, γ, μ, ω; sigdigits = 3, <keyword arguments>)

Find the smallest immune boosting parameter at which the endemic equilibrium switches 
    from stable to unstable.

`sigdigits` is the number of significant figures calculated. All other keyword arguments 
    are passed to `maxrealequileigen`.
"""
function findpsi(β, γ, μ, ω; sigdigits = 3, kwargs...)
    psivector = TRIALPSIS # see consts.jl
    return _findpsi(psivector, β, γ, μ, ω, sigdigits; kwargs...)
end 

function _findpsi(psivector, β, γ, μ, ω, sigdigits; kwargs...)
    # first iteraction 
    ind = findpsiindex(psivector, β, γ, μ, ω; kwargs...) 
    return _findpsi(ind, psivector, β, γ, μ, ω, sigdigits; kwargs...)
end 

_findpsi(ind::Nothing, psivector, β, γ, μ, ω, sigdigits; kwargs...) = NaN

function _findpsi(ind, psivector, β, γ, μ, ω, sigdigits; kwargs...)
    # vector for second iteration 
    v = _findpsivector(ind, psivector) 
    # 2nd to penultimate iterations
    for _ ∈ 1:(sigdigits - 2) 
        v = findpsivector(v, β, γ, μ, ω) 
    end 
    # final iteration
    ind = findpsiindex(v, β, γ, μ, ω)
    return v[ind]
end 

function findpsiindex(psivector, β, γ, μ, ω; kwargs...)
    eigs = [ realmaxequileigen(β, γ, μ, ψ, ω; kwargs...) for ψ ∈ psivector ]
    return findfirst(x -> x > 0, eigs)
end

function findpsivector(psivector, β, γ, μ, ω)
    ind = findpsiindex(psivector, β, γ, μ, ω) 
    return _findpsivector(ind, psivector)
end

function _findpsivector(ind, psivector)
    if ind == 1 lower = 0. 
    else        lower = psivector[ind-1]
    end
    upper = psivector[ind]
    diff = upper - lower 
    return collect(lower:diff/10:upper)
end
#=
function eigentracks(β, γ, μ, psis, ω)
    eigenvalues1 = equileigen(SirnsParameters(β, γ, μ, psis[1], ω))
    eigenvaluesmatrix = Matrix{ComplexF64}(undef, length(eigenvalues1), length(psis))
    for (i, ψ) ∈ enumerate(psis) 
        ev = equileigen(SirnsParameters(β, γ, μ, ψ, ω), sortby = λ -> (imag(λ),real(λ)))
        eigenvaluesmatrix[:, i] = ev
    end 
    return @ntuple eigenvalues1 eigenvaluesmatrix
end
=#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run simulations to find limit cycles for bifurcation plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
    bifurcationlimits(R0, γ, μ, psis, ω; <keyword arguments>)

Run simulations for a range of values of the immune boosting coefficient, and report 
    minimum and maximum proportions (after a warm-up period).

Outputs a `NamedTuple` containing `mins` and `maxs`
"""
function bifurcationlimits(R0, γ, μ, psis, ω; kwargs...)
    β = R0 * (γ + μ)
    mins = Vector{Float64}(undef, length(psis))
    maxs = Vector{Float64}(undef, length(psis))
    _bifurcationlimits!(mins, maxs, β, γ, μ, psis, ω; kwargs...)
    return @ntuple mins maxs
end

function _bifurcationlimits!(mins, maxs, β, γ, μ, psis, ω; kwargs...)
    for (i, ψ) ∈ enumerate(psis) 
        @unpack mn, mx = bifurcationlimit(β, γ, μ, ψ, ω; kwargs...)
        mins[i] = mn
        maxs[i] = mx
    end 
end

"""
    function bifurcationlimit(β, γ, μ, ψ, ω; <keyword arguments>)

Runs a simulation and reports the maximum and minimum proportions infectious.

Outputs a `NamedTuple` containing `mn`, the minimum value, and `mx`, the maximum value.

Default `tspan` for simulation is `( -1000., 10. )`. Maximum and minimum values are 
    calculated from the period with positive time (e.g. with the default, simulation
    runs for 1010 years and reports minimum and maximum values from final 10 years).

== Keyword arguments ==
* `maxiters = 1e5`: passed to `DifferentialEquations.solve` 
* `tspan = ( -1000., 10. )`: simulation duration 
* `S0 = .7`: Simulation initial proportion susceptible
* `I0 = .1`: Simulation initial proportion infectious
"""
function bifurcationlimit(β, γ, μ, ψ, ω; 
        maxiters = 1e5, tspan = ( -1000., 10. ), S0 = .7, I0 = .1
    )
    p = SirnsParameters(β, γ, μ, ψ, ω)
    u0 = sirns_u0(S0, I0; p, equalrs = true)
    sol = run_sirns(u0, p, tspan; maxiters)
    I = modelcompartments(sol, 2)
    mn = minimum(I)
    mx = maximum(I) 
    return @ntuple mn mx
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions for DrWatson.produce_or_load 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
    pl_bifurcationlimits(config)

Runs function `bifurcationlimits` in `DrWatson.produce_or_load`

`config` is a `NamedTuple` that contains parameters `R0`, `γ`, `μ`, `psis`, `ω`, 
    and optionally `maxiters`.
"""
function pl_bifurcationlimits(config)
    if haskey(config, :maxiters)
        @unpack R0, γ, μ, psis, ω, maxiters = config
    else 
        @unpack R0, γ, μ, psis, ω = config
        maxiters = 1e5 
    end
    minmax = bifurcationlimits(R0, γ, μ, psis, ω; maxiters) # returns an ntuple
    return ntuple2dict(minmax)
end
