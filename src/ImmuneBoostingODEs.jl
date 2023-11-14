
module ImmuneBoostingODEs 

using CairoMakie, CSV, DataFrames, DifferentialEquations, DrWatson, FFTW, Optimization, 
    SciMLSensitivity, StatsBase, Zygote
using LinearAlgebra: eigen
using OptimizationFlux: ADAM
using Roots: ZeroProblem
import Base: oneto

include("structs.jl")
include("consts.jl")
include("processdata.jl")
include("equilibria.jl")
include("model.jl")
include("analysedata.jl")
include("optimization.jl")
include("equilibriumsurface.jl")
include("plotformating.jl")

export 
    # structs.jl 
    SirnsParameters, LambdaParms,
    # consts.jl"
    COLOURVECTOR, COLOUR_I, COLOUR_R, COLOUR_S, MONTHDAYS,
    #processdata.jl
    printrawdate, processagedata, processcrgtvdata, processrsvdata,
    # equilibria.jl
    bifurcationlimits, equil, equileigen, equili, equilplotdata, equilr, 
    equilri, equils, findpsi, pl_bifurcationlimits, realmaxequileigen,
    # model.jl
    casespertimeblock, modelcompartments, pl_modelincidence, reducetransmission!, 
    restoretransmission!, run_sirns, sirns!, sirns_u0, 
    # analysedata.jl
    fourierhmdata,
    # optimization.jl 
    iterativeopt, opt_changetransmission!, opt_restoretransmission!, 
    # equilibriumsurface.jl
    labelequilibriumsurface!, plotequilibriumsurface!,
    # plotformating.jl 
    formataxis!, labelplots!, setvalue!

# Temporary code, awaiting update to Makie.jl 

import Makie: create_linepoints 
include("plottingfunction.jl")
export create_linepoints

# When Makie.jl is updated, Makie does not need to be an explicit dependency of this 
# project, and the file "plottingfunction.jl" can be deleted
 
end # module ImmuneBoostingODEs 
