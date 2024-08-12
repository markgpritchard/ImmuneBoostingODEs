
# Runs code from other scripts and produces figures for manuscript 

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
include("displayrsvanalysis.jl") # provides data, rsvsim_psi0, rsvsim_psi5, rsvsim_psi13_2,
    # agedata, rsvparms

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model flowchart / figure 1 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

introplot = Figure(; figure_padding = 1, size = ( 400, 200 ))
modelflowchart!(introplot)
# view the figure 
introplot
# save the figure
safesave(plotsdir("introplot.pdf"), introplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Effect of boosting on equilibria / figure 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

equilibriumplot = let 
    fig = Figure(; size=( 800, 750 ))
    plotequilibrium!(
        fig, critvector, bifurcationI_1_5, bifurcationI_5, 
        bifurcationI_10, bifurcationI_15, basic5, basic1_5, basic5_, equilparms
    )
    gd = GridLayout(fig[3, 1:3])
    plotsi!(
        gd, stablesim, 
        unstablein12, unstableout12, unstablein5, unstableout5, 
        limitcycle12, limitcycle5
    )
    rowsize!(fig.layout, 3, Auto(0.75))
    labelplots!("D", gd)
    fig
end

equilibriumplot
safesave(plotsdir("equilibriumplot.pdf"), equilibriumplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fourier transforms / figure 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fourierplot = Figure(; size = ( 800, 300 ))
plotfourier!(fourierplot, unforced6mfreqs, unforced400dfreqs, unforced25freqs, 
    unforced6mdensities, unforced400ddensities, unforced25densities, forced6mfreqs, 
    forced400dfreqs, forced25freqs, forced6mdensities, forced400ddensities, 
    forced25densities, simparms)

fourierplot
safesave(plotsdir("fourierplot.pdf"), fourierplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Respiratory syncytial virus data with simulations / figure 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fittedparametersfig = let 
    fig = Figure(; size=( 800, 1000 ))
    plotvvector = [ 
        plotvals01, plotvals02, plotvals04, plotvals1, plotvals2, plotvals4, plotvals10 
    ]
    parametervector = [
        rsvparameters01, rsvparameters02, rsvparameters04, 
        rsvparameters1, rsvparameters2, rsvparameters4, rsvparameters10 
    ]
    plotfittedsimulations!(fig, plotvvector, parametervector, data, crgtdata, saveat)  
    
    fig
end

safesave(plotsdir("fittedparametersfig.pdf"), fittedparametersfig)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stringency Index / supplementary figure 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stringencyplot = Figure(; size = ( 400, 400 ))
plotstringency!(stringencyplot, crgtdata, reduceday, increaseday) 

stringencyplot
safesave(plotsdir("stringencyplot.pdf"), stringencyplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proportions in each compartment at equilibrium / supplementary figure 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

equilibriumproportionsplot = Figure(; size = ( 400, 800 ))
plotequilibriumproportions!(equilibriumproportionsplot, equilSs, equilIs, equilRs, equilparms)

equilibriumproportionsplot
safesave(plotsdir("equilibriumproportionsplot.pdf"), equilibriumproportionsplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Non-pharmaceutical interventions plot / supplementary figure 3 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

npisimulationplot = Figure(; size = ( 800, 500 ))
plotnpi!(npisimulationplot, npisim_phi0, npisim_phi5, npisim_phi13_2, npiparms) 

npisimulationplot
safesave(plotsdir("npisimulationplot.pdf"), npisimulationplot)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-stratified respiratory syncytial virus data  / supplementary figure 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rsvagefigure = Figure(; size = ( 400, 600 ))
plotrsvage!(rsvagefigure, agedata)

rsvagefigure
safesave(plotsdir("rsvagefigure.pdf"), rsvagefigure)
