
using DrWatson
#using Revise
@quickactivate :ImmuneBoostingODEs
using CairoMakie, DifferentialEquations
CairoMakie.activate!() # allows figures to be saved as vector files


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
    ψ = 0

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, 1., 1.) 
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
    β1 = ϕ = 0
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot the Strigency Index 

stringencyplot = Figure(resolution = ( 400, 400 ))
let 
    ax = Axis(stringencyplot[1, 1])
    vspan!(ax, reduceday, increaseday, color = (:gray, 0.1))
    lines!(ax, crgtdata.Date, crgtdata.StringencyIndex_Average; color = COLOURVECTOR[1])
    Label(stringencyplot[2, 1], "Date"; fontsize = 11.84, tellwidth = false)
    Label(stringencyplot[1, 0], "Stringency index"; 
        fontsize = 11.84, rotation = π/2, tellheight = false)
    colgap!(stringencyplot.layout, 1, 5)
    rowgap!(stringencyplot.layout, 1, 5)
    formataxis!(ax)
end
stringencyplot
safesave(plotsdir("stringencyplot.pdf"), stringencyplot)

# Plot simulations 

npisimulationplot = Figure(resolution = ( 800, 500 ))

let 
    @unpack reductiontime = npiparms

    ga = GridLayout(npisimulationplot[1, 1])
    axs = [ Axis(ga[i, j]) for i ∈ [ 1, 3, 5 ], j ∈ [ 1, 3 ] ]
    for (i, m) ∈ enumerate([ npisim_phi0, npisim_phi5, npisim_phi13_2 ])
        @unpack cases, compartments = m 
        lines!(axs[i, 1], compartments[:gt][2:end], cases; color = COLOURVECTOR[3])
        lines!(axs[i, 2], compartments[:gt], compartments[:S]; color = COLOUR_S, label = "S")
        lines!(axs[i, 2], compartments[:gt], compartments[:R1]; color = COLOURVECTOR[4], label = "R₁")
        lines!(axs[i, 2], compartments[:gt], compartments[:R2]; color = COLOURVECTOR[5], label = "R₂")
        lines!(axs[i, 2], compartments[:gt], compartments[:R3]; color = COLOURVECTOR[6], label = "R₃")
        for j ∈ 1:2 vspan!(axs[i, j], reductiontime, reductiontime + 1, color = (:gray, 0.1)) end
   
    end
    leg = Legend(ga[0, 3], axs[2, 2])#; padding = ( 5, 5, 3, 3 ))

    linkxaxes!(axs...)
    for j ∈ 1:2 linkyaxes!(axs[:, j]...) end
    for i ∈ 1:3, j ∈ 1:2
        formataxis!(axs[i, j]; hidex = i != 3, hidexticks = i != 3)
        i != 3 && hidespines!(axs[i, j], :b)
    end
    formataxis!(leg)

    Label(ga[0, 1], "R₀ = 1.215 ± 10%, ψ = 0"; fontsize = 11.84, halign = :left, tellwidth = false)
    Label(ga[2, 1], "R₀ = 1.285 ± 8.2%, ψ = 5"; fontsize = 11.84, halign = :left, tellwidth = false)
    Label(ga[4, 1], "R₀ = 1.6, ψ = 13.2"; fontsize = 11.84, halign = :left, tellwidth = false)

    Label(ga[6, 1:3], "Time, years"; fontsize = 11.84, tellwidth = false)
    Label(ga[1:5, 0], "Weekly recorded incidence"; fontsize = 11.84, rotation = π/2, tellheight = false)
    Label(ga[1:5, 2], "Proportions"; fontsize = 11.84, rotation = π/2, tellheight = false)
    for r ∈ [ 1, 3, 5, 6 ] rowgap!(ga, r, 5) end
    for c ∈ [ 1, 3 ] colgap!(ga, c, 5) end
end 
npisimulationplot

safesave(plotsdir("npisimulationplot.pdf"), npisimulationplot)
