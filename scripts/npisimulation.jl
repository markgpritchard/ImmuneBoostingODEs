
using DrWatson
@quickactivate :ImmuneBoostingODEs
using DifferentialEquations


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple models without having a large 
# number of values in global scope potentially leading to unexpected results
npiparms = let 
    γ = 48.7    # generation time 7.5 days
    μ = .0087   # Scotland's birth rate = 48000 / 5.5e6
    tspan = ( -1000.65, 10. ) 
    βreduction = .8
    reductiontime = 5.2
    θ = .01 # proportion of incident infections recorded
    simulateddates = collect(-7/365:7/365:10)
    @ntuple γ μ tspan βreduction reductiontime θ simulateddates
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data from Oxford Covid-19 Government Response Tracker
crgtdata = processcrgtvdata("OxCGRT_compact_subnational_v1.csv", "crgt.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stringency Index in Scotland
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Use the Strigency Index to determine when the infection parameter is expected to change

# Find dates (as fractions of year) when Strigency Index goes above then below 50 
reduceday, increaseday = let 
    inds = findall(x -> x >= 50, crgtdata.StringencyIndex_Average)
    reduceday = crgtdata.Date[inds[1]]
    increaseday = crgtdata.Date[last(inds)]
    reduceday, increaseday
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define the ODE problem
prob = let 
    # different parameters may be used in the simulations
    @unpack γ, μ, tspan = npiparms 
    R0 = 1.6 
    β0 = R0 * (γ + μ)
    ω = 365.25 / 400
    β1 = .1 
    ϕ = -.25
    ψ = 0.0

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, β0, β0, β0) 
    I0 = .007
    S0 = .5
    u0 = sirns_u0(S0, I0; equalrs = true, p, t0 = -1000.65)
    ODEProblem(sirns!, u0, tspan, p)
end

## Run the simulations 

# Without immune boosting
npisim_phi0 = let 
    @unpack γ, μ, tspan, βreduction, reductiontime, θ, simulateddates = npiparms
    R0 = 1.215
    ψ = 0
    immuneduration = .5
    β1 = .1 
    ϕ = -.25
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    
    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, 1.) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = -.65)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reductiontime, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(reductiontime + 1, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)
    
    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = simulateddates, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 5e6)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

npisim_phi5 = let 
    @unpack γ, μ, tspan, βreduction, reductiontime, θ, simulateddates = npiparms
    R0 = 1.285
    ψ = 5
    immuneduration = .5
    β1 = .082 
    ϕ = -.25
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, 1.) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = -.65)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reductiontime, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(reductiontime + 1, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = simulateddates, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 5e6)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

npisim_phi13_2 = let 
    @unpack γ, μ, tspan, βreduction, reductiontime, θ, simulateddates = npiparms
    R0 = 1.6
    ψ = 13.2
    immuneduration = .5
    β1 = ϕ = 0.0
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, 1.)
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = -.65)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reductiontime, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(reductiontime + 1, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = simulateddates, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 5e6)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 
