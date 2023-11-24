
# Runs code from other scripts and produces figures formatted for poster presentation 

using DrWatson
@quickactivate :ImmuneBoostingODEs
using CairoMakie
CairoMakie.activate!() # allows figures to be saved as vector files

include("immuneduration.jl") # provides immunedurationmodel
include("equilibriumvalues.jl") # provides equilSs, equilIs, equilRs, equilparms, 
    # critvector, bifurcationI_1_5, bifurcationI_5, bifurcationI_10, bifurcationI_15, 
    # basic5, basic1_5, basic5_, stablesim, unstablein12, unstableout12, unstablein5, 
    # unstableout5, limitcycle12, limitcycle5
include("fouriertransforms.jl") # provides simparms, unforced6mfreqs, unforced400dfreqs, 
    # unforced25freqs, unforced6mdensities, unforced400ddensities, unforced25densities, 
    # forced6mfreqs, forced400dfreqs, forced25freqs
include("npisimulation.jl") # provides crgtdata, reduceday, increaseday, npisim_phi0, 
    # npisim_phi5, npisim_phi13_2, npiparms


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model flowchart
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

introplot = Figure(; figure_padding = 1, size = ( 400, 200 ))
modelflowchart!(introplot)
# view the figure 
introplot
# save the figure
safesave(plotsdir("posterintroplot.pdf"), introplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proportions in infectious compartment at equilibrium
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

equilibriumproportionsplot = Figure(; figure_padding = ( 10, 50, 10, 0 ), size = ( 300, 250 ))
plotequilibriumproportions!(equilibriumproportionsplot, equilIs, equilparms)

equilibriumproportionsplot
safesave(plotsdir("posterequilibriumproportionsplot.pdf"), equilibriumproportionsplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Effect of boosting on equilibria
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

equilibriumplot = Figure(; size = ( 600, 400 ))

plotequilibriuma!(equilibriumplot, critvector, equilparms; labelplot = false)
plotequilibriumb!(equilibriumplot, critvector, bifurcationI_1_5, bifurcationI_5, 
    bifurcationI_10, bifurcationI_15, equilparms; labelplot = false)
let 
    gc = GridLayout(equilibriumplot[2, 1:3])
    plotsi!(gc, stablesim, unstablein12, unstableout12, unstablein5, unstableout5, 
        limitcycle12, limitcycle5
    )
end
rowsize!(equilibriumplot.layout, 1, Auto(1.3))

equilibriumplot
safesave(plotsdir("posterequilibriumplot.pdf"), equilibriumplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation without seasonal forcing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unforcedplot = Figure(; size = ( 600, 225 ))
let 
    gl = GridLayout(unforcedplot[1, 1])
    plotequilibriumc!(gl, basic5, basic1_5, basic5_)
end

unforcedplot
safesave(plotsdir("posterunforcedplot.pdf"), unforcedplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Non-pharmaceutical interventions plot / Simulation with seasonal forcing 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

npisimulationplot = Figure(; size = ( 600, 375 ))
plotnpi!(npisimulationplot, npisim_phi0, npisim_phi5, npisim_phi13_2, npiparms) 

npisimulationplot
safesave(plotsdir("posternpisimulationplot.pdf"), npisimulationplot)
