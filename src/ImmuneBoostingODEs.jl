
module ImmuneBoostingODEs 

import CSV
import NaNMath

using CairoMakie: @L_str, Axis, Axis3, Figure, GridLayout, Label, RGBf, lines!, surface!
using DataFrames: DataFrame, insertcols!, leftjoin!, rename!, select!, subset, subset!
using Dates: day, month, year
using DifferentialEquations: ODEProblem, Vern9, solve
using DrWatson: @dict, @ntuple, @unpack
using DrWatson: datadir, load, ntuple2dict, produce_or_load, tostringdict
using FFTW: fft
using LineSearches: Static
using LinearAlgebra: eigen
using Optim: LBFGS
using Roots: ZeroProblem
using SciMLBase: successful_retcode
using StatsBase: mean, quantile, std
using Turing: @addlogprob!, @model
using Turing: Beta, Exponential, InitFromParams, InitFromPrior, LogNormal, Normal, Uniform
using Turing: arraydist, maximum_a_posteriori, maximum_likelihood, pdf, truncated

include("consts.jl")
include("processdata.jl")
include("equilibria.jl")
include("model.jl")
include("rsvfitmodel.jl")
include("analysedata.jl")
include("equilibriumsurface.jl")
include("plotting.jl")

## structs.jl 
export SirnsParameters, LambdaParms
## consts.jl"
export COLOURVECTOR, COLOUR_I, COLOUR_R, COLOUR_S
## processdata.jl
export printrawdate, processagedata, processcrgtvdata, processmobilitydata, processrsvdata
## equilibria.jl
export bifurcationlimits 
export equil, equileigen, equili, equilplotdata, equilr, equilri, equils 
export findpsi 
export pl_bifurcationlimits 
export realmaxequileigen
## model.jl
export casespertimeblock 
export modelcompartments 
export pl_modelincidence 
export reducetransmission! 
export restoretransmission! 
export run_sirns
export sirns!
export sirns_u0
## rsvfitmodel.jl
export fitmodel, initialparams_map, initialparams_mle
## analysedata.jl
export fittedsimulationquantiles 
export fittedsimulationsetup
export fourierhmdata, simulatedfourierhmdata
export loadrsvdata
export memosolver
export runfittedsimulations
## equilibriumsurface.jl
export labelequilibriumsurface!, plotequilibriumsurface!
## plotting.jl
export plotchains
export plotequilibriumc!
export plotfittedsimulationquantiles!
export plotfourier!
export plotnpi!
export plotrsvage!
export plotsi!
export plotstringency!
 
end  # module ImmuneBoostingODEs 
