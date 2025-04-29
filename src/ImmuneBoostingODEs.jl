
module ImmuneBoostingODEs 

using CairoMakie, CSV, DataFrames, DifferentialEquations, DrWatson, FFTW, ForwardDiff
using PlotFormatting, StatsBase
using LinearAlgebra: eigen
using Memoization: @memoize
using Roots: ZeroProblem
import Base: ==, hash

include("structs.jl")
include("consts.jl")
include("processdata.jl")
include("equilibria.jl")
include("model.jl")
include("analysedata.jl")
include("equilibriumsurface.jl")
include("plotting.jl")

export 
    ## structs.jl 
    SirnsParameters, 
    LambdaParms,
    ## consts.jl"
    COLOURVECTOR, 
    COLOUR_I, 
    COLOUR_R, 
    COLOUR_S, 
    MONTHDAYS,
    ## processdata.jl
    printrawdate, 
    processagedata, 
    processcrgtvdata, 
    processrsvdata,
    ## equilibria.jl
    bifurcationlimits, 
    equil, 
    equileigen, 
    equili, 
    equilplotdata, 
    equilr, 
    equilri, 
    equils, 
    findpsi, 
    pl_bifurcationlimits, 
    realmaxequileigen,
    ## model.jl
    casespertimeblock, 
    modelcompartments, 
    pl_modelincidence, 
    reducetransmission!, 
    restoretransmission!, 
    run_sirns, 
    sirns!, 
    sirns_u0, 
    ## analysedata.jl
    fittedsimulationquantiles, 
    fittedsimulationsetup, 
    fourierhmdata, 
    loadrsvdata, 
    memosolver, 
    runfittedsimulations, 
    ## equilibriumsurface.jl
    labelequilibriumsurface!, 
    plotequilibriumsurface!,
    ## plotting.jl
    plotchains, 
    plotequilibriumc!, 
    plotfittedsimulationquantiles!, 
    plotfourier!, 
    plotnpi!, 
    plotrsvage!, 
    plotsi!, 
    plotstringency!
 
end  # module ImmuneBoostingODEs 
