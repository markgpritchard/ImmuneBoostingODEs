
module ImmuneBoostingODEs 

using CairoMakie, CSV, DataFrames, DifferentialEquations, DrWatson, FFTW, StatsBase
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
include("plotformating.jl")

export 
    # structs.jl 
    SirnsParameters, LambdaParms,
    # consts.jl"
    COLOURVECTOR, COLOUR_I, COLOUR_R, COLOUR_S, MONTHDAYS,
    #processdata.jl
    printrawdate, processagedata, processcrgtvdata, processrsvdata,
    # equilibria.jl
    bifurcationlimits, equil, equileigen, equili, equilplotdata, equilr, equilri, equils, 
    findpsi, pl_bifurcationlimits, realmaxequileigen,
    # model.jl
    casespertimeblock, modelcompartments, pl_modelincidence, reducetransmission!, 
    restoretransmission!, run_sirns, sirns!, sirns_u0, 
    # analysedata.jl
    fourierhmdata, memosolver,
    # equilibriumsurface.jl
    labelequilibriumsurface!, plotequilibriumsurface!,
    # plotting.jl
    modelflowchart!, plotchains, plotequilibrium!, plotequilibriuma!, plotequilibriumb!, 
    plotequilibriumc!,plotequilibriumproportions!, plotfourier!, plotimmuneduration!, 
    plotnpi!, plotrsvage!, plotrsvsim!, plotsi!, plotstringency!,
    # plotformating.jl 
    formataxis!, labelplots!, setvalue!
 
end  # module ImmuneBoostingODEs 
