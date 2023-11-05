
module ImmuneBoostingODEs 

using CairoMakie, CSV, DataFrames, DifferentialEquations, DrWatson, FFTW, StatsBase
using LinearAlgebra: eigen
using Roots: ZeroProblem
import Base: oneto

include("structs.jl")
include("consts.jl")
include("processdata.jl")
include("equilibria.jl")
include("model.jl")
include("analysedata.jl")
include("equilibriumsurface.jl")
include("plotformating.jl")

export 
    # structs.jl 
    SirnsParameters, LambdaParms,
    # consts.jl"
    COLOURVECTOR, COLOUR_I, COLOUR_R, COLOUR_S, MONTHDAYS,
    #processdata.jl
    processagedata, processcrgtvdata, processrsvdata,
    # equilibria.jl
    bifurcationlimits, equil, equileigen, equili, equilplotdata, equilr, 
    equilri, equils, findpsi, pl_bifurcationlimits, realmaxequileigen,
    # model.jl
    casespertimeblock, modelcompartments, pl_modelincidence, reducetransmission!, 
    restoretransmission!, run_sirns, sirns!, sirns_u0, 
    # analysedata.jl
    fourierhmdata,
    # equilibriumsurface.jl
    labelequilibriumsurface!, plotequilibriumsurface!,
    # plotformating.jl 
    formataxis!, labelplots!, setvalue!

end # module ImmuneBoostingODEs 
