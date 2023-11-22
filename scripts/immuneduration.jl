
# Duration of immunity without immune boosting

using DrWatson
@quickactivate :ImmuneBoostingODEs
using CairoMakie
CairoMakie.activate!() # allows figures to be saved as vector files

# Generate duration of immunity by assuming waning at rate 1 and no new infections
immunedurationplot = let 
    p = SirnsParameters(.0, .0, .0, .0, .0, .0, 1.)
    u0 = sirns_u0(.0, .0; equalrs = false, p)
    sol = run_sirns(u0, p, ( 0., 2.5 ))
    mc = modelcompartments(sol, p)

    immunedurationplot = Figure(; size = ( 400, 350 ))
    axs = [ Axis(immunedurationplot[i, 1]) for i ∈ 1:2 ]
    lines!(axs[1], mc[:gt], mc[:Rtotal]; color = COLOUR_R)
    # rate of leaving immunity is 3ω * R3
    lines!(axs[2], mc[:gt], 3 .* mc[:R3]; color = COLOUR_R)
    for i ∈ 1:2 vlines!(axs[i], 1; color = :black, linestyle = :dot) end
    Label(immunedurationplot[1, 0], "Proportion immune"; fontsize = 11.84, rotation = π/2, tellheight = false)
    Label(immunedurationplot[2, 0], "Rate of return\nto susceptibility"; fontsize = 11.84, rotation = π/2, tellheight = false)
    Label(immunedurationplot[3, 1], "Time, multiples of mean immune duration"; fontsize = 11.84, tellwidth = false)
    formataxis!(axs[1]; hidespines = ( :r, :t, :b ), hidex = true, hidexticks = true)
    formataxis!(axs[2])
    colgap!(immunedurationplot.layout, 1, 5)
    rowgap!(immunedurationplot.layout, 2, 5)

    immunedurationplot
end

safesave(plotsdir("immunedurationplot.pdf"), immunedurationplot)
