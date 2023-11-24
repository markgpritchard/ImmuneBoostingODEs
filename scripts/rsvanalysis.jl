
using DrWatson
@quickactivate :ImmuneBoostingODEs
using DataFrames, DifferentialEquations

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RSV data from Scotland
data = processrsvdata("respiratory_scot.csv", "rsv.csv")

# Age-specific data
agedata = processagedata("respiratory_age.csv", "rsv_age.csv")

# Data from Oxford Covid-19 Government Response Tracker
crgtdata = processcrgtvdata("OxCGRT_compact_subnational_v1.csv", "crgt.csv")

# To avoid splitting outbreaks, count cases from April each year 
let 
    april1value = MONTHDAYS[4] / 365
    offsetdate = data.Date .- april1value
    aprilyear = round.(Int, offsetdate, RoundDown)
    aprilfractiondate = offsetdate - aprilyear
    insertcols!(data, :AprilYear => aprilyear)
    insertcols!(data, :AprilFractionDate => aprilfractiondate)
    # insert cumulative cases since last April 
    cumulativecases = Vector{Float64}(undef, size(data, 1))
    cumulativecases[1] = data.Cases[1]
    for i ∈ axes(data, 1)
        i == 1 && continue
        if data.AprilYear[i] == data.AprilYear[i-1]
            cumulativecases[i] = data.Cases[i] + cumulativecases[i-1]
        else 
            cumulativecases[i] = data.Cases[i]
        end 
    end 
    insertcols!(data, :AprilCumulativeCases => cumulativecases)
end 
println("In the 12 months from 1 April each year")
for y ∈ 2016:2022 
    inds = findall(x -> y <= x < y + 1, data.AprilYear)
    println("    $(sum(data.Cases[inds])) cases in $y")
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple models without having a large 
# number of values in global scope potentially leading to unexpected results
rsvparms = let 
    γ = 48.7    # generation time 7.5 days
    μ = .0087   # Scotland's birth rate = 48000 / 5.5e6
    ϕ = -.5π

    # When is the infection parameter expected to change?
    # Find dates (as fractions of year) when Strigency Index goes above then below 50
    inds = findall(x -> x >= 50, crgtdata.StringencyIndex_Average)
    reduceday = crgtdata.Date[inds[1]]
    increaseday = crgtdata.Date[last(inds)]
    println("Stringency ≥ 50 on $(printrawdate(crgtdata.RawDate[inds[1]]))")
    println("Stringency < 50 on $(printrawdate(crgtdata.RawDate[last(inds)]))")
    # note the Stringency Index is plotted by code in `npisimulation.jl`

    # Vector of recorded cases 
    casesvector = data.Cases 

    ## Times when we have data pre-lockdown (i.e. times to save simulation)
    # add one pre-data date so that we can calculate a weekly incidence for the first data point
    savetimes = [ minimum(data.Date) - 7 / 365; data.Date ] # 354 elements

    @ntuple ϕ γ μ casesvector reduceday increaseday savetimes
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Optimize magnitude of the effect of non-pharmaceutical interventions 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Callbacks for optimization 
opt_cbs = let 
    @unpack increaseday, reduceday = rsvparms
    save_positions = ( false, false )
    resetcb = PresetTimeCallback(reduceday + 1e-9, opt_restoretransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, opt_changetransmission!; save_positions)
    CallbackSet(resetcb, rescb)
end

multipliers_psi0 = let 
    @unpack casesvector, increaseday, reduceday, savetimes = rsvparms
    R0 = 1.215
    β1 = .1
    ψ = 0
    initialvalues = [ .8, 1. ]
    iterativeopt(R0, β1, ψ, initialvalues; 
        casesvector, increaseday, opt_cbs, reduceday, savetimes) 
end

multipliers_psi5 = let 
    @unpack casesvector, increaseday, reduceday, savetimes = rsvparms
    R0 = 1.285
    β1 = .082
    ψ = 5
    initialvalues = [ .8, 1. ]
    iterativeopt(R0, β1, ψ, initialvalues;
        casesvector, increaseday, opt_cbs, reduceday, savetimes) 
end

multipliers_psi13_2 = let 
    @unpack casesvector, increaseday, reduceday, savetimes = rsvparms
    R0 = 1.6
    β1 = 0
    ψ = 13.2
    initialvalues = [ .8, 1. ]
    iterativeopt(R0, β1, ψ, initialvalues; 
        casesvector, increaseday, opt_cbs, reduceday, savetimes) 
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define the ODE problem 
prob = let 
    # different parameters may be used in the simulations
    @unpack ϕ, γ, μ = rsvparms
    R0 = 1.6 
    β0 = R0 * (γ + μ)
    ω = 365.25 / 400
    β1 = .0 
    ψ = 0

    tspan = ( 1015.35, 2023.6 )
    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, 1., 1.) 
    I0 = .007
    S0 = .5
    u0 = sirns_u0(S0, I0; equalrs = true, p, t0 = 1015.35)
    ODEProblem(sirns!, u0, tspan, p)
end

## Run the simulations 

# Without immune boosting
rsvsim_psi0 = let 
    @unpack ϕ, γ, μ, reduceday, increaseday, savetimes = rsvparms
    βreduction, βreturn = multipliers_psi0 
    R0 = 1.215
    ψ = 0
    immuneduration = .5
    β1 = .1
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    
    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = .35)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)
    
    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 1e7)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

rsvsim_psi5 = let 
    @unpack ϕ, γ, μ, reduceday, increaseday, savetimes = rsvparms
    βreduction, βreturn = multipliers_psi5
    R0 = 1.285
    ψ = 5
    immuneduration = .5
    β1 = .082
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = .35) 
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 1e7)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

rsvsim_psi13_2 = let 
    @unpack ϕ, γ, μ, reduceday, increaseday, savetimes = rsvparms
    βreduction, βreturn = multipliers_psi13_2
    R0 = 1.6
    ψ = 13.2
    immuneduration = .5
    β1 = ϕ = 0
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = .35)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 5e6)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 
