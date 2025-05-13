
module ImmuneBoostingODEs 

import Base: ==, hash

using AutoHashEquals: @auto_hash_equals
using CSV
using CairoMakie
using DataFrames
using Dates
using DifferentialEquations
using DrWatson
using FFTW
using ForwardDiff
using LinearAlgebra: eigen
using Memoization: @memoize
using PlotFormatting
using RollingFunctions: rollmean
using Roots: ZeroProblem
using StatsBase

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
    processmobilitydata,
    processrsvdata,
    ## equilibria.jl
    bifurcationlimits, 
    equil, 
    equileigen, 
    equili, 
    #equilplotdata, 
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
    sirns_u0_transformedp,
    transformparameters,
    transformedsirns!,
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
