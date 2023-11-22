
# Find and plot equilibrium values for model 

using DrWatson
@quickactivate :ImmuneBoostingODEs
using CairoMakie
CairoMakie.activate!() # allows figures to be saved as vector files

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple simulations, and across running 
# and plotting simulations, without having a large number of values in global scope 
# potentially leading to unexpected or inconsistent results

equilparms = let 
    γ       = 48.7    # generation time 7.5 days
    μ       = .0087   # Scotland's birth rate = 48_000 / 5.5e6
    psis    = collect(0:.1:20)
    ω       = .913    # mean duration of immunity without boosting 400 days 
    durations = [ collect(.01:.01:2.5); collect(2.5:.1:80) ] # Used for plotting 
    # effect of changing mean duration of immunity. Set to provide finer resolution 
    # for first 2.5 years
    R0s     = collect(.0125:.025:15)
    @ntuple γ μ psis ω durations R0s
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Equilibrium proportion in each compartment for different parameter values 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Values of I^* for different values of R0 and ϕ, which is plotted in Figure 1

let 
    @unpack γ, μ, psis, ω, R0s = equilparms
    betas = R0s .* (γ + μ)
    
    global equilSs = zeros(length(betas), length(psis))
    global equilIs = zeros(length(betas), length(psis))
    global equilRs = zeros(length(betas), length(psis))
    for (i, β) ∈ enumerate(betas), (j, ϕ) ∈ enumerate(psis)
        I = equili(SirnsParameters(β, γ, μ, ϕ, ω))
        if 0 < I <= 1   # endemic equilibrium
            equilSs[i, j] = equils(SirnsParameters(β, γ, μ, ϕ, ω))
            equilIs[i, j] = I
            equilRs[i, j] = equilr(SirnsParameters(β, γ, μ, ϕ, ω), I) 
        else            # disease-free equilibium
            equilSs[i, j] = 1.
            equilIs[i, j] = .0
            equilRs[i, j] = .0
        end 
    end 
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stability of endemic equilibria for different values of ω 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For Figure 2 we want the stability for a range of values of immune duration at 
# different R0 values

critvector = let 
    @unpack γ, μ, psis, durations = equilparms
    R0s = [ 1.6, 5, 10, 15 ]
    betas = R0s .* (γ + μ)
    omegas = 1 ./ durations
    [ [ realmaxequileigen(SirnsParameters(β, γ, μ, ψ, ω); warntol = Inf) 
            for ψ ∈ psis, ω ∈ omegas ]
        for β ∈ betas ]
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations to show stability and limit cycles
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Examples with R0 = 1.6 are for Figure 1, and examples with other R0 values for 
# supplementary figures.
 
# find the relevant critical value of ϕ 
let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    β = R0 * (γ + μ)
    println("ϕ = $(findpsi(β, γ, μ, ω; warntol = Inf))")
end     # ϕ = 1.18

# stable simulation
stablesim = let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    ϕ = .5
    β = R0 * (γ + μ)
    p = SirnsParameters(β, γ, μ, ϕ, ω)
    u0 = sirns_u0(.8, .03; p, equalrs = true)
    tspan = ( 0., 35. ) # it has to run a long time to avoid looking like a limit cycle
    sol = run_sirns(u0, p, tspan)
    modelcompartments(sol, p)
end

# unstable simulations
let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    ϕ = 1.5
    β = R0 * (γ + μ)
    p = SirnsParameters(β, γ, μ, ϕ, ω)
    
    # find the limit cycle 
    let
        u0 = sirns_u0(.7, .1; p, equalrs = true)
        tspan = ( -1000., 10. )
        sol = run_sirns(u0, p, tspan)
        global limitcycle12 = modelcompartments(sol, p)
    end

    # approach the cycle from inside
    let
        u0 = sirns_u0(.64, .008; p, equalrs = true)
        tspan = ( 0., 3.8 )
        sol = run_sirns(u0, p, tspan)
        global unstablein12 = modelcompartments(sol, p)
    end

    # approach the cycle from outside
    let
        u0 = sirns_u0(.8, .03; p, equalrs = true)
        tspan = ( 0., 3.1 )
        sol = run_sirns(u0, p, tspan)
        global unstableout12 = modelcompartments(sol, p)
    end 
end 

let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    ϕ = 5
    β = R0 * (γ + μ)
    p = SirnsParameters(β, γ, μ, ϕ, ω)
    
    # find the limit cycle 
    let
        u0 = sirns_u0(.7, .1; p, equalrs = true)
        tspan = ( -1000., 10. )
        sol = run_sirns(u0, p, tspan)
        global limitcycle5 = modelcompartments(sol, p)
    end

    # approach the cycle from inside
    let
        u0 = sirns_u0(.75, .01; p, equalrs = true)
        tspan = ( 0., 3.2 )
        sol = run_sirns(u0, p, tspan)
        global unstablein5 = modelcompartments(sol, p)
    end

    # approach the cycle from outside
    let
        u0 = sirns_u0(.8, .03; p, equalrs = true)
        tspan = ( 0., 3.5 )
        sol = run_sirns(u0, p, tspan)
        global unstableout5 = modelcompartments(sol, p)
    end 
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

let 
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    β = R0 * (γ + μ)
    tspan = ( 0., 10. )

    let 
        ϕ = .5 
        p = SirnsParameters(β, γ, μ, ϕ, ω)
        u0 = sirns_u0(.7, .01; p, equalrs = true)
        sol = run_sirns(u0, p, tspan)
        global basic5 = modelcompartments(sol, p)
    end 

    let 
        ϕ = 1.5 
        p = SirnsParameters(β, γ, μ, ϕ, ω)
        u0 = sirns_u0(.7, .01; p, equalrs = true)
        sol = run_sirns(u0, p, tspan)
        global basic1_5 = modelcompartments(sol, p)
    end 

    let 
        ϕ = 5 
        p = SirnsParameters(β, γ, μ, ϕ, ω)
        u0 = sirns_u0(.7, .01; p, equalrs = true)
        sol = run_sirns(u0, p, tspan)
        global basic5_ = modelcompartments(sol, p)
    end 
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Values for bifurcation plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bifurcationI_1_5 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 1.5
    config = @dict R0 γ μ psis ω
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end

bifurcationI_5 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 5
    config = @dict R0 γ μ psis ω maxiters = 1e6
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end

bifurcationI_10 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 10
    config = @dict R0 γ μ psis ω maxiters = 1e6
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end

bifurcationI_15 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 15
    config = @dict R0 γ μ psis ω maxiters = 1e6
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots of equilibrium values 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Proportions in each compartment at equilibrium 

equilibriumproportionsplot = Figure(; size = ( 400, 800 ))
let 
    @unpack γ, μ, psis, ω, R0s = equilparms
    betas = R0s .* (γ + μ)

    axs = [ Axis3(equilibriumproportionsplot[i, 1]; azimuth = .7π) for i ∈ 1:3 ]
    values = [ equilSs, equilIs, equilRs ]
    keywords = [ "susceptible", "infectious", "recovered" ]
    for (ax, val, kw) ∈ zip(axs, values, keywords) 
        plotequilibriumsurface!(ax, R0s, psis, val; rasterize = 30)
        labelequilibriumsurface!(ax, [ 
            "Basic reproduction\nnumber, R₀", 
            "Boosting parameter, ψ", 
            "Equilibrium\nproportion $kw" 
        ];
        xlaboffset = 25, ylaboffset = 25, zlaboffset = 80)
    end
    labelplots!([ "A", "B", "C" ], equilibriumproportionsplot; cols = 1, rows = [ 1, 2, 3 ]) 
end 

equilibriumproportionsplot
# Save the plot
safesave(plotsdir("equilibriumproportionsplot.pdf"), equilibriumproportionsplot)

## Threshold values for stable equilibrium 

equilibriumplot = Figure(; size = ( 800, 550 ))

let 
    @unpack psis, durations = equilparms
    R0s = [ 1.5, 5, 10, 15 ]
    R0lbls = [ "1.5", "5", "10", "15" ]
    colourvec = [ i == 1 ? :black : ( COLOURVECTOR[i+2], .5 ) for i ∈ eachindex(critvector) ]
    ga = GridLayout(equilibriumplot[1, 1])

    # Main plot (<= 2.5 years of immunity)
    xinds = findall(x -> x <= 2.5, durations)
    yinds = findall(x -> x <= 10, psis)

    mainax = Axis(ga[1, 1])
    for (i, fullcrit) ∈ enumerate(critvector)
        crit = fullcrit[yinds, xinds]
        contour!(mainax, durations[xinds], psis[yinds], crit'; 
            color = colourvec[i], levels = [ 0 ])
    end 
    text!(mainax, 1.5, 3; text = "Unstable", fontsize = 11.84, align = ( :center, :center ))
    text!(mainax, .4, .5; text = "Stable", fontsize = 11.84, align = ( :center, :center ))

    # Inset plot for longer duration immunity 
    miniax = Axis(ga[1, 1]; 
        halign = :right, valign = :top, height = Relative(.5), width = Relative(.37))
    for (i, crit) ∈ enumerate(critvector)
        contour!(miniax, durations, psis, crit'; color = colourvec[i], levels = [ 0 ])
    end 

    formataxis!(mainax; setorigin = true)
    formataxis!(miniax; setorigin = true, hidespines = nothing, trimspines = false)
    Label(ga[2, 1], "Mean duration of immunity, years"; fontsize = 11.84, tellwidth = false)
    Label(ga[1, 0], "Boosting coefficient, ψ*"; 
        fontsize = 11.84, rotation = π/2, tellheight = false)
    colgap!(ga, 1, 5)
    rowgap!(ga, 1, 5)
    labelplots!("A", ga; rows = 1) 

    leg = Legend(equilibriumplot[1, 3], [ LineElement(color = c) for c ∈ colourvec ], R0lbls, "R₀:"; 
        titlefont = :regular, orientation = :vertical)
    formataxis!(leg; horizontal = false)
end

# Bifurcation diagram
let 
    colourvec = [ i == 1 ? :black : ( COLOURVECTOR[i+2], .5 ) for i ∈ eachindex(critvector) ]
    @unpack psis = equilparms 

    gb = GridLayout(equilibriumplot[1, 2])
    axs = [ Axis(gb[i, 1]; yticks = LinearTicks(3)) for i ∈ 1:4 ]
    for (i, bf) ∈ enumerate([ bifurcationI_1_5, bifurcationI_5, bifurcationI_10, bifurcationI_15])
        @unpack mins, maxs = bf
        lines!(axs[i], psis, mins; color = colourvec[i])
        lines!(axs[i], psis, maxs; color = colourvec[i])
        formataxis!(axs[i]; hidex = i != 4, hidexticks = i != 4)
        if i != 4 hidespines!(axs[i], :b) end
    end 

    Label(gb[5, 1], "Boosting coefficient, ψ"; fontsize = 11.84, tellwidth = false)
    Label(gb[1:4, 0], "Proportion infectious"; fontsize = 11.84, rotation = π/2, tellheight = false)

    colgap!(gb, 1, 5)
    rowgap!(gb, 4, 5)
    labelplots!("B", gb; rows = 1) 
end 

# Compartments over time 
let 
    gc = GridLayout(equilibriumplot[2, 1:3])

    Label(gc[0, 1], "ψ = 0.5"; halign = :left, fontsize = 11.84, tellwidth = false)
    axs1 = [ Axis(gc[i, 1]) for i ∈ 1:3 ]
    lines!(axs1[1], basic5[:gt], basic5[:S]; color = COLOUR_S)
    lines!(axs1[2], basic5[:gt], basic5[:I]; color = COLOUR_I)
    lines!(axs1[3], basic5[:gt], basic5[:R1]; color = COLOURVECTOR[4])
    lines!(axs1[3], basic5[:gt], basic5[:R2]; color = COLOURVECTOR[5])
    lines!(axs1[3], basic5[:gt], basic5[:R3]; color = COLOURVECTOR[6])
    
    Label(gc[0, 2], "ψ = 1.5"; halign = :left, fontsize = 11.84, tellwidth = false)
    axs2 = [ Axis(gc[i, 2]) for i ∈ 1:3 ]
    lines!(axs2[1], basic1_5[:gt], basic1_5[:S]; color = COLOUR_S)
    lines!(axs2[2], basic1_5[:gt], basic1_5[:I]; color = COLOUR_I)
    lines!(axs2[3], basic1_5[:gt], basic1_5[:R1]; color = COLOURVECTOR[4])
    lines!(axs2[3], basic1_5[:gt], basic1_5[:R2]; color = COLOURVECTOR[5])
    lines!(axs2[3], basic1_5[:gt], basic1_5[:R3]; color = COLOURVECTOR[6])
    
    Label(gc[0, 3], "ψ = 5"; halign = :left, fontsize = 11.84, tellwidth = false)
    axs3 = [ Axis(gc[i, 3]) for i ∈ 1:3 ]
    lines!(axs3[1], basic5_[:gt], basic5_[:S]; color = COLOUR_S)
    lines!(axs3[2], basic5_[:gt], basic5_[:I]; color = COLOUR_I)
    lines!(axs3[3], basic5_[:gt], basic5_[:R1]; color = COLOURVECTOR[4])
    lines!(axs3[3], basic5_[:gt], basic5_[:R2]; color = COLOURVECTOR[5])
    lines!(axs3[3], basic5_[:gt], basic5_[:R3]; color = COLOURVECTOR[6])
    
    leg = Legend(gc[3, 4], [ LineElement(color = c) for c ∈ COLOURVECTOR[4:6] ], 
        [ "$ℓ" for ℓ ∈ [ "R₁", "R₂", "R₃" ]]; 
        titlefont = :regular, orientation = :vertical)

    linkxaxes!(axs1..., axs2..., axs3...)
    for j ∈ 1:3 
        linkyaxes!(axs1[j], axs2[j], axs3[j])
        for (i, axs) ∈ enumerate([ axs1, axs2, axs3 ])
            formataxis!(axs[j]; 
                hidex = j != 3, hidexticks = j != 3, hidey = i != 1, hideyticks = i != 1)
            if i != 1 hidespines!(axs[j], :l) end
            if j != 3 hidespines!(axs[j], :b) end
        end 
    end 
    formataxis!(leg; horizontal = false)

    Label(gc[4, 1:3], "Time (years)"; fontsize = 11.84, tellwidth = false)
    Label(gc[1, 0], "Proportion\nsusceptible"; fontsize = 11.84, tellheight = false)
    Label(gc[2, 0], "Proportion\ninfectious"; fontsize = 11.84, tellheight = false)
    Label(gc[3, 0], "Proportions\nresistant"; fontsize = 11.84, tellheight = false)
    
    labelplots!("C", gc; rows = 1) 

    for c ∈ [ 1, 4 ] colgap!(gc, c, 5) end
    for r ∈ [ 1, 4 ] rowgap!(gc, r, 5) end
    for r ∈ [ 2, 3 ] rowgap!(gc, r, 7) end
end 
colsize!(equilibriumplot.layout, 1, Auto(1.2))
colgap!(equilibriumplot.layout, 2, 5)

equilibriumplot
# Save the plot
safesave(plotsdir("equilibriumplot.pdf"), equilibriumplot)

## Susceptible--infectious planes

siplot = Figure(; size = ( 800, 250 ))
let
    ax1 = Axis(siplot[1, 1])
    lines!(ax1, stablesim[:S], stablesim[:I]; color = COLOURVECTOR[4])
    for ℓ ∈ [ 3050, 5900 ]
        x = stablesim[:S][ℓ]
        y = stablesim[:I][ℓ]
        u = x - stablesim[:S][ℓ-1] 
        v = y - stablesim[:I][ℓ-1] 
        arrows!(ax1, [x], [y], [u], [v]; color = COLOURVECTOR[4])
    end

    ax2 = Axis(siplot[1, 2])
    lines!(ax2, limitcycle12[:S], limitcycle12[:I]; color = :black)
    for (i, r) ∈ enumerate([ unstablein12, unstableout12 ])
        lines!(ax2, r[:S], r[:I]; color = COLOURVECTOR[i+4])
        ℓ = length(r[:S])
        x = r[:S][ℓ]
        y = r[:I][ℓ]
        u = x - r[:S][ℓ-1] 
        v = y - r[:I][ℓ-1] 
        arrows!(ax2, [x], [y], [u], [v]; color = COLOURVECTOR[i+4])
    end

    ax3 = Axis(siplot[1, 3])
    lines!(ax3, limitcycle5[:S], limitcycle5[:I]; color = :black)
    for (i, r) ∈ enumerate([ unstablein5, unstableout5 ])
        lines!(ax3, r[:S], r[:I]; color = COLOURVECTOR[i+4])
        ℓ = length(r[:S])
        x = r[:S][ℓ]
        y = r[:I][ℓ]
        u = x - r[:S][ℓ-1] 
        v = y - r[:I][ℓ-1] 
        arrows!(ax3, [x], [y], [u], [v]; color = COLOURVECTOR[i+4])
    end

    linkaxes!(ax1, ax2, ax3)

    formataxis!(ax1)
    formataxis!(ax2; hidey = true, hideyticks = true, hidespines = (:l, :r, :t))
    formataxis!(ax3; hidey = true, hideyticks = true, hidespines = (:l, :r, :t))
end
Label(siplot[0, 1], "ψ = 0.5"; halign = :left, fontsize = 11.84, tellwidth = false)
Label(siplot[0, 2], "ψ = 1.5"; halign = :left, fontsize = 11.84, tellwidth = false)
Label(siplot[0, 3], "ψ = 5"; halign = :left, fontsize = 11.84, tellwidth = false)
Label(siplot[2, 1:3], "Proportion susceptible"; fontsize = 11.84, tellwidth = false)
Label(siplot[1, 0], "Proportion infectious"; fontsize = 11.84, rotation = π/2, tellheight = false)
colgap!(siplot.layout, 1, 5) 
for r ∈ 1:2 rowgap!(siplot.layout, r, 5) end

siplot
# Save the plot
safesave(plotsdir("siplot.pdf"), siplot)
