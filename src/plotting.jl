
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model flowchart  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function modelflowchart!(fig)
    fontsize = 11.84 # font size for parameter labels on arrows 

    ga = GridLayout(fig[1, 1])
    ax = Axis(ga[1, 1]) # axis that figure will be plotted on

    for i ∈ 1:5
        color = ( [ COLOUR_S, COLOUR_I, COLOUR_R, COLOUR_R, COLOUR_R ][i], .2 )
        label = [ "S", "I", "R₁", "R₂", "R₃" ][i]
        # generate coloured boxes for each compartment
        poly!(ax, Rect(i - .25, -.25, .5, .5); color)
        # label each box
        text!(ax, i, 0; text = label, align = ( :center, :center ))
    end 
    # arrows between boxes
    for x ∈ 2:4 arrows!(ax, [ x + .28 ], [ 0 ], [ .4 ], [ 0 ]; color = :black) end
    arrows!(ax, [ .52 ], [ .52 ], [ .2 ], [ -.2 ]; color = :black)
    for x ∈ 1:5 arrows!(ax, [ x + .28 ], [ -.28 ], [ .2 ], [ -.2 ]; color = :black) end
    # dashed lines
    xs1 = collect(1.28:.05:1.63); ys1 = zeros(length(xs1))
    linesegments!(ax, xs1, ys1; color = :red)
    arrows!(ax, [ 1.63 ], [ 0 ], [ .05 ], [ 0 ]; color = :red)
    xs2 = collect(4:-.05:3.05); ys2 = @. .52 - (xs2 - 3.5)^2
    linesegments!(ax, xs2, ys2; color = :red)
    xs3 = collect(5:-.05:3.05); ys3 = @. .57 - .3 * (xs3 - 4)^2
    linesegments!(ax, xs3, ys3; color = :red)
    arrows!(ax, [ 3.05 ], [ .2975 ], [ -.05 ], [ -.0275 ]; color = :red)
    lines!(ax, [ 5, 5 ], [ -.28, -.93 ]; color = :black)
    lines!(ax, [ 1, 5 ], [ -.93, -.93 ]; color = :black)
    arrows!(ax, [ 1 ], [ -.93 ], [ 0 ], [ .58 ]; color = :black)
    # label the arrows / lines
    text!(ax, 1.5, .1; text = "β(t) I(t)", color = :red, fontsize, align = ( :center, :bottom ))
    text!(ax, 2.5, .1; text = "γ", color = :black, fontsize, align = ( :center, :bottom ))
    for x ∈ 3:4
        text!(ax, x + .5, .1; 
            text = "3ω", color = :black, fontsize, align = ( :center, :bottom ))
    end
    text!(ax, .72, .45; text = "μ", color = :black, fontsize, align = ( :left, :bottom ))
    for x ∈ 1:5
        text!(ax, x + .48, -.35; text = "μ", color = :black, fontsize, align = ( :left, :bottom ))
    end
    text!(ax, 3.8, .65; text = "ψ β(t) I(t)", color = :red, fontsize, align = ( :center, :bottom ))
    text!(ax, 3, -1.03; text = "3ω", color = :black, fontsize, align = ( :center, :top ))

    formataxis!(ax; hidex = true, hidexticks = true, hidey = true, hideyticks = true, 
        hidespines = ( :l, :t, :r, :b ))
    setvalue!(ax, (.5, -1.1 )); setvalue!(ax, (5.5 , .75 ))
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Effect of boosting on equilibria 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotequilibrium!(fig, critvector, bifurcationI_1_5, bifurcationI_5, 
        bifurcationI_10, bifurcationI_15, basic5, basic1_5, basic5_, equilparms
    )
    plotequilibriuma!(fig, critvector, equilparms)
    plotequilibriumb!(fig, critvector, bifurcationI_1_5, bifurcationI_5, bifurcationI_10, 
        bifurcationI_15, equilparms)
    plotequilibriumc!(fig, basic5, basic1_5, basic5_)
    colsize!(fig.layout, 1, Auto(1.2))
    colgap!(fig.layout, 2, 5)
end

function plotequilibriuma!(fig, critvector, equilparms; labelplot = true)
    @unpack psis, durations = equilparms
    R0s = [ 1.5, 5, 10, 15 ]
    R0lbls = [ "1.5", "5", "10", "15" ]
    colourvec = [ i == 1 ? :black : ( COLOURVECTOR[i+2], .5 ) for i ∈ eachindex(critvector) ]
    ga = GridLayout(fig[1, 1])

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
    labelplot && labelplots!("A", ga; rows = 1) 

    leg = Legend(fig[1, 3], [ LineElement(color = c) for c ∈ colourvec ], R0lbls, "R₀:"; 
        titlefont = :regular, orientation = :vertical)
    formataxis!(leg; horizontal = false)
end

# Bifurcation diagram

function plotequilibriumb!(fig, critvector, bifurcationI_1_5, bifurcationI_5, bifurcationI_10, 
        bifurcationI_15, equilparms; labelplot = true
    )
    colourvec = [ i == 1 ? :black : ( COLOURVECTOR[i+2], .5 ) for i ∈ eachindex(critvector) ]
    @unpack psis = equilparms 

    gb = GridLayout(fig[1, 2])
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
    labelplot && labelplots!("B", gb; rows = 1) 
end 

# Compartments over time 

function plotequilibriumc!(fig::Figure, basic5, basic1_5, basic5_)
    gc = GridLayout(fig[2, 1:3])
    plotequilibriumc!(gc, basic5, basic1_5, basic5_)
    labelplots!("C", gc; rows = 1) 
end 

function plotequilibriumc!(gc::GridLayout, basic5, basic1_5, basic5_)
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
    
    for c ∈ [ 1, 4 ] colgap!(gc, c, 5) end
    for r ∈ [ 1, 4 ] rowgap!(gc, r, 5) end
    for r ∈ [ 2, 3 ] rowgap!(gc, r, 7) end
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proportions in each compartment at equilibrium
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotequilibriumproportions!(fig, equilSs, equilIs, equilRs, equilparms)
    axs = [ Axis3(fig[i, 1]; azimuth = .7π) for i ∈ 1:3 ]
    values = [ equilSs, equilIs, equilRs ]
    keywords = [ "susceptible", "infectious", "recovered" ]
    _plotequilibriumproportions!(fig, axs, values, keywords, equilparms)
    labelplots!([ "A", "B", "C" ], fig; cols = 1, rows = [ 1, 2, 3 ]) 
end 

function plotequilibriumproportions!(fig, equilIs, equilparms)
    axs = [ Axis3(fig[i, 1]; azimuth = .7π) for i = 1 ]
    values = [ equilIs ]
    keywords = [ "infectious" ]
    _plotequilibriumproportions!(fig, axs, values, keywords, equilparms)
end 

function _plotequilibriumproportions!(fig, axs, values, keywords, equilparms)
    @unpack γ, μ, psis, ω, R0s = equilparms
    betas = R0s .* (γ + μ)


    for (ax, val, kw) ∈ zip(axs, values, keywords) 
        plotequilibriumsurface!(ax, R0s, psis, val; rasterize = 30)
        labelequilibriumsurface!(ax, [ 
            "Basic reproduction\nnumber, R₀", 
            "Boosting parameter, ψ", 
            "Equilibrium\nproportion $kw" 
        ];
        xlaboffset = 25, ylaboffset = 25, zlaboffset = 80)
    end
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fourier transforms 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotfourier!(fig, unforced6mfreqs, unforced400dfreqs, unforced25freqs, 
    unforced6mdensities, unforced400ddensities, unforced25densities, forced6mfreqs, 
    forced400dfreqs, forced25freqs, forced6mdensities, forced400ddensities, 
    forced25densities, simparms
)
gl = GridLayout(fig[1, 1])

Label(gl[0, 1], "Immune duration 6 months"; fontsize = 11.84, halign = :left, tellwidth = false)
Label(gl[0, 2], "1.1 years"; fontsize = 11.84, halign = :left, tellwidth = false)
Label(gl[0, 3], "2.5 years"; fontsize = 11.84, halign = :left, tellwidth = false)

_plotfouriera!(gl, unforced6mfreqs, unforced400dfreqs, unforced25freqs, unforced6mdensities, 
    unforced400ddensities, unforced25densities, simparms)
_plotfourierb!(gl, forced6mfreqs, forced400dfreqs, forced25freqs, forced6mdensities, 
    forced400ddensities, forced25densities, simparms) 

labelplots!([ "A", "B"], gl; rows = [ 0, 2])
for r ∈ [ 1, 3 ] rowgap!(gl, r, 5) end
for c ∈ [ 1, 5 ] colgap!(gl, c, 5) end
end

function _plotfouriera!(gl, unforced6mfreqs, unforced400dfreqs, unforced25freqs, 
    unforced6mdensities, unforced400ddensities, unforced25densities, simparms
)
@unpack psis = simparms
freqvector = [ unforced6mfreqs, unforced400dfreqs, unforced25freqs ]
densityvector = [ unforced6mdensities, unforced400ddensities, unforced25densities ]

axs = [ Axis(gl[1, i]; yscale = log) for i ∈ 1:3 ]
i = 1
for (ax, freqs, densities) ∈ zip(axs, freqvector, densityvector)
    inds = findall(x -> .1 <= x <= 365 / 21, freqs) 
    hm = heatmap!(ax, psis, freqs[inds], log10.(densities[inds, :]'); 
        colormap = :CMRmap, colorrange = ( -6, 0 ))
    ax.yticks = ([.2, .5, 1, 2, 4, 365 / 35], [ "5 y", "2 y", "1 y", "6 m", "3 m", "5 w"])
    formataxis!(ax; hidex = true, hidey = i != 1)
    if i == 1 
        co = Colorbar(gl[1:2, 4], hm) 
        formataxis!(co)
    end
    i += 1
end
Label(gl[1:2, 0], "Period"; 
    fontsize = 11.84, rotation = π/2, tellheight = false)
Label(gl[1:2, 5], "log₁₀ spectral density"; 
    fontsize = 11.84, rotation = -π/2, tellheight = false)
end

function _plotfourierb!(gl, forced6mfreqs, forced400dfreqs, forced25freqs, forced6mdensities, 
    forced400ddensities, forced25densities, simparms
) 
@unpack psis = simparms
freqvector = [ forced6mfreqs, forced400dfreqs, forced25freqs ]
densityvector = [ forced6mdensities, forced400ddensities, forced25densities ]

axs = [ Axis(gl[2, i]; yscale = log) for i ∈ 1:3 ]
i = 1
for (ax, freqs, densities) ∈ zip(axs, freqvector, densityvector)
    inds = findall(x -> .1 <= x <= 365 / 21, freqs) 
    hm = heatmap!(ax, psis, freqs[inds], log10.(densities[inds, :]'); 
        colormap = :CMRmap, colorrange = ( -6, 0 ))
    ax.yticks = ([.2, .5, 1, 2, 4, 365 / 35], [ "5 y", "2 y", "1 y", "6 m", "3 m", "5 w"])
    formataxis!(ax; hidey = i != 1)
    i += 1
end
Label(gl[3, 1:3], "Boosting coefficient, ψ"; 
    fontsize = 11.84, tellwidth = false)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Immune duration plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotimmuneduration!(fig, modeloutput)
    axs = [ Axis(fig[i, 1]) for i ∈ 1:2 ]
    lines!(axs[1], modeloutput[:gt], modeloutput[:Rtotal]; color = COLOUR_R)
    # rate of leaving immunity is 3ω * R3
    lines!(axs[2], modeloutput[:gt], 3 .* modeloutput[:R3]; color = COLOUR_R)
    for i ∈ 1:2 vlines!(axs[i], 1; color = :black, linestyle = :dot) end
    Label(fig[1, 0], "Proportion immune"; fontsize = 11.84, rotation = π/2, tellheight = false)
    Label(fig[2, 0], "Rate of return\nto susceptibility"; fontsize = 11.84, rotation = π/2, tellheight = false)
    Label(fig[3, 1], "Time, multiples of mean immune duration"; fontsize = 11.84, tellwidth = false)
    formataxis!(axs[1]; hidespines = ( :r, :t, :b ), hidex = true, hidexticks = true)
    formataxis!(axs[2])
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Non-pharmaceutical interventions plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotnpi!(fig, npisim_phi0, npisim_phi5, npisim_phi13_2, npiparms) 
    @unpack reductiontime = npiparms

    ga = GridLayout(fig[1, 1])
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

    Label(
        ga[0, 1], "R₀ = 1.215 with seasonal forcing ± 10%, ψ = 0"; 
        fontsize = 11.84, halign = :left, tellwidth = false
    )
    Label(
        ga[2, 1], "R₀ = 1.285 with seasonal forcing ± 8.2%, ψ = 5"; 
        fontsize = 11.84, halign = :left, tellwidth = false
    )
    Label(
        ga[4, 1], "R₀ = 1.6 with no seasonal forcing, ψ = 13.2"; 
        fontsize = 11.84, halign = :left, tellwidth = false
    )

    Label(ga[6, 1:3], "Time, years"; fontsize = 11.84, tellwidth = false)
    Label(ga[1:5, 0], "Weekly incidence"; fontsize = 11.84, rotation = π/2, tellheight = false)
    Label(ga[1:5, 2], "Proportions"; fontsize = 11.84, rotation = π/2, tellheight = false)
    for r ∈ [ 1, 3, 5, 6 ] rowgap!(ga, r, 5) end
    for c ∈ [ 1, 3 ] colgap!(ga, c, 5) end
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-stratified respiratory syncytial virus data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotrsvage!(fig, agedata)
    colours = [ ( :gray, .4); [ COLOURVECTOR[i] for i ∈ [ 5, 6, 4 ] ] ]

    axs = [ Axis(fig[i, 1]) for i ∈ 1:7 ]
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
    leg = Legend(fig[0, 1], lineelements, labels, "Year:"; titlefont = :regular)
    formataxis!(leg)

    Label(fig[8, 1], "Month"; fontsize = 11.84, tellwidth = false)
    Label(fig[1:7, 0], "Annual cumulative incidence"; fontsize = 11.84, rotation = π/2, tellheight = false)
    for r ∈ [ 1, 8 ] rowgap!(fig.layout, r, 5) end
    for c ∈ [ 1 ] colgap!(fig.layout, c, 5) end
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Respiratory syncytial virus data with simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotrsvsim!(fig, data, rsvsim_psi0, rsvsim_psi5, rsvsim_psi13_2, rsvparms) 
    @unpack reduceday, increaseday = rsvparms
    ga = GridLayout(fig[1, 1])
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Susceptible--infectious planes 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotsi!(fig::Figure, stablesim, unstablein12, unstableout12, unstablein5, 
        unstableout5, limitcycle12, limitcycle5
    )
    _plotsi!(fig, stablesim, unstablein12, unstableout12, unstablein5, unstableout5, 
        limitcycle12, limitcycle5)
    colgap!(fig.layout, 1, 5) 
    for r ∈ 1:2 rowgap!(fig.layout, r, 5) end
end

function plotsi!(fig::GridLayout, stablesim, unstablein12, unstableout12, unstablein5, 
        unstableout5, limitcycle12, limitcycle5
    )
    _plotsi!(fig, stablesim, unstablein12, unstableout12, unstablein5, unstableout5, 
        limitcycle12, limitcycle5)
    colgap!(fig, 1, 5) 
    for r ∈ 1:2 rowgap!(fig, r, 5) end
end

function _plotsi!(fig, stablesim, unstablein12, unstableout12, unstablein5, 
        unstableout5, limitcycle12, limitcycle5
    )
    ax1 = Axis(fig[1, 1])
    lines!(ax1, stablesim[:S], stablesim[:I]; color = COLOURVECTOR[4])
    for ℓ ∈ [ 3050, 5900 ]
        x = stablesim[:S][ℓ]
        y = stablesim[:I][ℓ]
        u = x - stablesim[:S][ℓ-1] 
        v = y - stablesim[:I][ℓ-1] 
        arrows!(ax1, [x], [y], [u], [v]; color = COLOURVECTOR[4])
    end

    ax2 = Axis(fig[1, 2])
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

    ax3 = Axis(fig[1, 3])
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
    Label(fig[0, 1], "ψ = 0.5"; halign = :left, fontsize = 11.84, tellwidth = false)
    Label(fig[0, 2], "ψ = 1.5"; halign = :left, fontsize = 11.84, tellwidth = false)
    Label(fig[0, 3], "ψ = 5"; halign = :left, fontsize = 11.84, tellwidth = false)
    Label(fig[2, 1:3], "Proportion susceptible"; fontsize = 11.84, tellwidth = false)
    Label(fig[1, 0], "Proportion infectious"; fontsize = 11.84, rotation = π/2, tellheight = false)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Strigency Index
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotstringency!(fig, crgtdata, reduceday, increaseday) 
    ax = Axis(fig[1, 1])
    vspan!(ax, reduceday, increaseday, color = (:gray, 0.1))
    lines!(ax, crgtdata.Date, crgtdata.StringencyIndex_Average; color = COLOURVECTOR[1])
    Label(fig[2, 1], "Date"; fontsize = 11.84, tellwidth = false)
    Label(fig[1, 0], "Stringency index"; 
        fontsize = 11.84, rotation = π/2, tellheight = false)
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    formataxis!(ax)
end
