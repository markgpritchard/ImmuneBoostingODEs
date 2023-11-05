
using DrWatson
#using Revise
@quickactivate :ImmuneBoostingODEs
using CairoMakie, DataFrames, DifferentialEquations
#import StatsPlots: plot as spplot # CairoMakie also exports a function called `plot`
CairoMakie.activate!() # allows figures to be saved as vector files


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

# When is the infection parameter expected to change?

# Find dates (as fractions of year) when Strigency Index goes above then below 40 
reduceday, increaseday = let 
    inds = findall(x -> x >= 40, crgtdata.StringencyIndex_Average)
    reduceday = crgtdata.Date[inds[1]]
    increaseday = crgtdata.Date[last(inds)]
    reduceday, increaseday
end

# note the Stringency Index is plotted by code in `npisimulation.jl`


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple models without having a large 
# number of values in global scope potentially leading to unexpected results
rsvparms = let 
    γ = 48.7    # generation time 7.5 days
    μ = .0087   # Scotland's birth rate = 48000 / 5.5e6
    tspan = ( 1015.35, 2023.6 ) 
    ϕ = -.5π
    n_adapts = 1000
    samples = 2000
    n_chains = 3
    acceptance = .65
    cases = data.Cases

    ## Times when we have data pre-lockdown (i.e. times to save simulation)
    # add one pre-data date so that we can calculate a weekly incidence for the first data point
    savetimes = [ minimum(data.Date) - 7 / 365; data.Date ] # 354 elements
    @ntuple ϕ γ μ acceptance cases n_adapts n_chains samples savetimes tspan
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define the ODE problem 
prob = let 
    # different parameters may be used in the simulations
    @unpack ϕ, γ, μ, tspan = rsvparms
    R0 = 1.6 
    β0 = R0 * (γ + μ)
    ω = 365.25 / 400
    β1 = .0 
    ψ = 0

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, 1., 1.) 
    I0 = .007
    S0 = .5
    u0 = sirns_u0(S0, I0; equalrs = true, p, t0 = 1015.35)
    ODEProblem(sirns!, u0, tspan, p)
end

## Run the simulations 

# Without immune boosting
rsvsim_psi0 = let 
    @unpack ϕ, γ, μ, savetimes, tspan = rsvparms
    R0 = 1.215
    ψ = 0
    immuneduration = .5
    β1 = .1
    βreduction = .87 
    βreturn = 1.
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    
    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = -.65)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday - 30 / 365, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)
    
    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 1e7)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

rsvsim_psi5 = let 
    @unpack ϕ, γ, μ, savetimes, tspan = rsvparms
    R0 = 1.285
    ψ = 5
    immuneduration = .5
    β1 = .082
    βreduction = .85
    βreturn = .95
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = -.65)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday - 30 / 365, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 1e7)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

rsvsim_psi13_2 = let 
    @unpack ϕ, γ, μ, savetimes, tspan = rsvparms
    R0 = 1.6
    ψ = 13.2
    immuneduration = .5
    β1 = ϕ = 0
    βreduction = .74
    βreturn = .825
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = -.65)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday - 30 / 365, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 5e6)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot data with simulations 

rsvfigure = Figure(resolution = ( 400, 600 ))

let 
    ga = GridLayout(rsvfigure[1, 1])
    ax1 = Axis(ga[1, 1])
    lines!(ax1, data.Date, data.Cases; color = :black)
    vspan!(ax1, reduceday, increaseday, color = ( :gray, 0.1 )) 

    axs = [ Axis(ga[i, 1]) for i ∈ [ 3, 5, 7 ] ]
    for (i, m) ∈ enumerate([ rsvsim_psi0, rsvsim_psi5, rsvsim_psi13_2 ])
        @unpack cases, compartments = m 
        scatter!(axs[i], data.Date, data.Cases; color = :black, markersize = 3)
        lines!(axs[i], compartments[:gt][2:end], cases; color = COLOURVECTOR[3])
        vspan!(axs[i], reduceday, increaseday, color = ( :gray, 0.1 )) 
    end

    linkaxes!(ax1, axs...)
    ax1.xticks = collect(2017:1:2023)
    formataxis!(ax1; hidex = true, hidexticks = true, hidespines = [ :b, :t, :r ])
    for i ∈ 1:3
        axs[i].xticks = collect(2017:1:2023)
        formataxis!(axs[i]; hidex = i != 3, hidexticks = i != 3)
        i != 3 && hidespines!(axs[i], :b)
    end

    Label(ga[0, 1], "Recorded cases in Scotland"; fontsize = 11.84, halign = :left, tellwidth = false)
    Label(ga[2, 1], "R₀ = 1.215 ± 10%, ψ = 0"; fontsize = 11.84, halign = :left, tellwidth = false)
    Label(ga[4, 1], "R₀ = 1.285 ± 8.2%, ψ = 5"; fontsize = 11.84, halign = :left, tellwidth = false)
    Label(ga[6, 1], "R₀ = 1.6, ψ = 13.2"; fontsize = 11.84, halign = :left, tellwidth = false)

    Label(ga[8, 1], "Time, years"; fontsize = 11.84, tellwidth = false)
    Label(ga[1:7, 0], "Weekly recorded incidence"; fontsize = 11.84, rotation = π/2, tellheight = false)
    for r ∈ [ 1, 3, 5, 7, 8 ] rowgap!(ga, r, 5) end
    for c ∈ [ 1 ] colgap!(ga, c, 5) end
end 
rsvfigure

safesave(plotsdir("rsvfigure.pdf"), rsvfigure)


## Plot age-stratified data

rsvagefigure = Figure(resolution = ( 400, 600 ))

let 
    colours = [ ( :gray, .4); [ COLOURVECTOR[i] for i ∈ [ 5, 6, 4 ] ] ]

    axs = [ Axis(rsvagefigure[i, 1]) for i ∈ 1:7 ]
    for (i, y) ∈ enumerate(collect(2016:1:2022))
        d = subset(agedata, :Year => x -> x.== y)
        for (j, a) ∈ enumerate(unique(agedata.AgeGroup))
            d2 = subset(d, :AgeGroup => x -> x.== a )
            if y <= 2019 color = colours[1] else color = colours[i-3] end
            lines!(axs[j], d2.FractionDate, getproperty(d2, Symbol("Rate$a")); color, label = "$y") 
            if y == 2022 
                maxy = maximum(getproperty(d2, Symbol("Rate$a")))
                text!(axs[j], 0, .9 * maxy; text = "$a", fontsize = 11.84, align = ( :left, :top ))
            end
        end
    end 

    for i ∈ 1:7
        formataxis!(axs[i]; hidex = i != 7, hidexticks = i != 7) 
        if i != 7 hidespines!(axs[i], :b) end
    end

    linkxaxes!(axs...)
    axs[7].xticks = (
        [ 0, .25, .5, .75, 1 ], 
        [ "April", "July", "Oct", "Jan", "April" ]
    )

    # create legend manually 
    lineelements = [ LineElement(; color = c) for c ∈ colours ]
    labels = [ "2016-19", "2020", "2021", "2022" ]
    leg = Legend(rsvagefigure[0, 1], lineelements, labels, "Year:"; titlefont = :regular)
    formataxis!(leg)

    Label(rsvagefigure[8, 1], "Month"; fontsize = 11.84, tellwidth = false)
    Label(rsvagefigure[1:7, 0], "Annual cumulative incidence"; fontsize = 11.84, rotation = π/2, tellheight = false)
    for r ∈ [ 1, 8 ] rowgap!(rsvagefigure.layout, r, 5) end
    for c ∈ [ 1 ] colgap!(rsvagefigure.layout, c, 5) end
end
rsvagefigure

safesave(plotsdir("rsvagefigure.pdf"), rsvagefigure)
