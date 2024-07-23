
using DrWatson
@quickactivate :ImmuneBoostingODEs
using DataFrames, DifferentialEquations

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RSV data from Scotland
data = processrsvdata("respiratory_scot.csv", "rsv.csv")

# Age-specific data
agedata = processagedata("respiratory_age.csv", "rsv_age.csv")

# Data from Oxford Covid-19 Government Response Tracker
crgtdata = processcrgtvdata("OxCGRT_compact_subnational_v1.csv", "crgt.csv")

# To avoid splitting outbreaks, count cases from April each year 
let 
    april1value = MONTHDAYS[4] / 365
    offsetdate = data.Date .- april1value
    aprilyear = round.(Int, offsetdate, RoundDown)
    aprilfractiondate = offsetdate - aprilyear
    insertcols!(data, :AprilYear => aprilyear)
    insertcols!(data, :AprilFractionDate => aprilfractiondate)
    # insert cumulative cases since last April 
    cumulativecases = Vector{Float64}(undef, size(data, 1))
    cumulativecases[1] = data.Cases[1]
    for i ∈ axes(data, 1)
        i == 1 && continue
        if data.AprilYear[i] == data.AprilYear[i-1]
            cumulativecases[i] = data.Cases[i] + cumulativecases[i-1]
        else 
            cumulativecases[i] = data.Cases[i]
        end 
    end 
    insertcols!(data, :AprilCumulativeCases => cumulativecases)
end 
println("In the 12 months from 1 April each year")
for y ∈ 2016:2022 
    inds = findall(x -> y <= x < y + 1, data.AprilYear)
    println("    $(sum(data.Cases[inds])) cases in $y")
end 

###

import ImmuneBoostingODEs: AbstractParameters 
using CSV
using CairoMakie
using Distributions
using Turing
using Pigeons
using Random

using Memoization


mutable struct SirrrsParameters{S} <: AbstractParameters where S
    τ   :: S
    c   :: Matrix{Float64}
    γ   :: Float64
    ν   :: Float64
    μ   :: Vector{Float64} 
    ψ   :: S
    ω   :: Float64
end     

function sirrrs!(du, u, p, t)
    # columns are age-groups, rows are compartments

    # if solely dividing by N_i, then when there is only one infant, if they happen to be 
    # infectious then they have an infeasible force of infection on the whole population. 
    # Therefore set a minimum N of one month's births. Will want to modify this slightly to 
    # make it continuously differentiable 
    #λ = [ 
    #    p.τ * (p.c[i, 1] * u[2, 1] / p.ν + sum([ p.c[i, j] * u[2, j] / sum(u[1:5, j]) for j ∈ 2:16 ])) * 365.25 
    #    for i ∈ 1:16
    #]
    λ = [ 
        p.τ * (
            p.c[i, 1] * sum(@view u[2, 1:5]) / (p.ν + sum(@view u[1:5, 1:5])) + 
            p.c[i, 2] * sum(@view u[2, 6:10]) / sum(@view u[1:5, 6:10]) + 
            p.c[i, 3] * sum(@view u[2, 11:15]) / sum(@view u[1:5, 11:15]) + 
            sum([ p.c[i, j] * u[2, (j + 12)] / sum(@view u[1:5, (j + 12)]) for j ∈ 4:16 ])
        ) * 365.25
        for i ∈ 1:16
    ]

    # age 0
    du[1, 1] = p.ν + 3 * p.ω * u[5, 1] - λ[1] * u[1, 1]                   # S
    du[2, 1] = λ[1] * u[1, 1] - p.γ * u[2, 1]                                 # I
    du[3, 1] = p.γ * u[2, 1] + λ[1] * p.ψ * (u[4, 1] + u[5, 1]) - (3 * p.ω) * u[3, 1]    # R1
    du[4, 1] = 3 * p.ω * u[3, 1] - (3 * p.ω + λ[1] * p.ψ) * u[4, 1]           # R2
    du[5, 1] = 3 * p.ω * u[4, 1] - (3 * p.ω + λ[1] * p.ψ) * u[5, 1]           # R3
    du[6, 1] = λ[1] * u[1, 1]  # cumulative cases

    for i ∈ 2:15  # ages 1 to 14 
        if i <= 5 
            j = 1 
        elseif i <= 10 
            j = 2 
        else 
            j = 3 
        end
        du[1, i] = 3 * p.ω * u[5, i] - λ[j] * u[1, i]                   # S
        du[2, i] = λ[j] * u[1, i] - p.γ * u[2, i]                                 # I
        du[3, i] = p.γ * u[2, i] + λ[j] * p.ψ * (u[4, i] + u[5, i]) - (3 * p.ω) * u[3, i]    # R1
        du[4, i] = 3 * p.ω * u[3, i] - (3 * p.ω + λ[j] * p.ψ) * u[4, i]           # R2
        du[5, i] = 3 * p.ω * u[4, i] - (3 * p.ω + λ[j] * p.ψ) * u[5, i]           # R3
        du[6, i] = λ[j] * u[1, i]  # cumulative cases
    end

    # ages 15 to 19 
    du[1, 16] = 3 * p.ω * u[5, 16] - λ[4] * u[1, 16] - (0.2 + p.μ[16]) * u[1, 16]                 # S
    du[2, 16] = λ[4] * u[1, 16] - p.γ * u[2, 16] - (0.2 + p.μ[16]) * u[2, 16]                                   # I
    du[3, 16] = p.γ * u[2, 16] + λ[4] * p.ψ * (u[4, 16] + u[5, 16]) - (3 * p.ω) * u[3, 16] - (0.2 + p.μ[16]) * u[3, 16]      # R1
    du[4, 16] = 3 * p.ω * u[3, 16] - (3 * p.ω + λ[4] * p.ψ) * u[4, 16] - (0.2 + p.μ[16]) * u[4, 16]             # R2
    du[5, 16] = 3 * p.ω * u[4, 16] - (3 * p.ω + λ[4] * p.ψ) * u[5, 16] - (0.2 + p.μ[16]) * u[5, 16]             # R3
    du[6, 16] = λ[4] * u[1, 16]  # cumulative cases

    for i ∈ 17:27  # ages 20 to 74 
        j = i - 12
        du[1, i] = 3 * p.ω * u[5, i] - λ[j] * u[1, i] + 0.2 * u[1, i-1] - (0.2 + p.μ[i]) * u[1, i]                 # S
        du[2, i] = λ[j] * u[1, i] - p.γ * u[2, i] + 0.2 * u[2, i-1] - (0.2 + p.μ[i]) * u[2, i]                                   # I
        du[3, i] = p.γ * u[2, i] + λ[j] * p.ψ * (u[4, i] + u[5, i]) - (3 * p.ω) * u[3, i] + 0.2 * u[3, i-1] - (0.2 + p.μ[i]) * u[3, i]      # R1
        du[4, i] = 3 * p.ω * u[3, i] - (3 * p.ω + λ[j] * p.ψ) * u[4, i] + 0.2 * u[4, i-1] - (0.2 + p.μ[i]) * u[4, i]             # R2
        du[5, i] = 3 * p.ω * u[4, i] - (3 * p.ω + λ[j] * p.ψ) * u[5, i] + 0.2 * u[5, i-1] - (0.2 + p.μ[i]) * u[5, i]             # R3
        du[6, i] = λ[j] * u[1, i]  # cumulative cases
    end

    # age >= 75 
    du[1, 28] = 3 * p.ω * u[5, 28] - λ[16] * u[1, 28] + 0.2 * u[1, 27] - p.μ[28] * u[1, 28]                 # S
    du[2, 28] = λ[16] * u[1, 28] - p.γ * u[2, 28] + 0.2 * u[2, 27] - p.μ[28] * u[2, 28]                                   # I
    du[3, 28] = p.γ * u[2, 28] + λ[16] * p.ψ * (u[4, 28] + u[5, 28]) - (3 * p.ω) * u[3, 28] + 0.2 * u[3, 27] - p.μ[28] * u[3, 28]      # R1
    du[4, 28] = 3 * p.ω * u[3, 28] - (3 * p.ω + λ[16] * p.ψ) * u[4, 28] + 0.2 * u[4, 27] - p.μ[28] * u[4, 28]             # R2
    du[5, 28] = 3 * p.ω * u[4, 28] - (3 * p.ω + λ[16] * p.ψ) * u[5, 28] + 0.2 * u[5, 27] - p.μ[28] * u[5, 28]             # R3
    du[6, 28] = λ[16] * u[1, 28]  # cumulative cases


end
#=
function _contactmatrixrowcolconversion(x)
    if x <= 15 
        return 1 + round(Int, (x - 1) / 5, RoundDown)
    else 
        return x - 12 
    end
end
=#
#=
_contactmatrixrowcolconversion(x) = x

_cmrc(x) = _contactmatrixrowcolconversion(x)
=#

function makecontactmatrix(source)
    contacts_df = CSV.read(source, DataFrame; header=false)
    #=contacts_raw = Matrix(contacts_df)
    contacts = Matrix{Float64}(undef, 28, 28)
    for i ∈ 1:16, j ∈ 1:16 
        contacts[i, j] = contacts_raw[_cmrc(i), _cmrc(j)]
    end
    return contacts =#
    return Matrix(contacts_df)
end

function makeoutputmatrices(sol)
    t = sol.t 
    susceptible = makeoutputmatrix_susceptible(sol)
    infectious = makeoutputmatrix_infectious(sol)
    resistant = makeoutputmatrix_resistant(sol)
    N = makeoutputmatrix_N(sol)
    cumcases = makeoutputmatrix_cumcases(sol)
    incidentcases = makeoutputmatrix_incidentcases(sol)
    return @ntuple t susceptible infectious resistant N cumcases incidentcases
end

function _makeoutputmatrix(sol::SciMLBase.AbstractODESolution{T, <:Any, <:Any}, ind) where T
    mat = Matrix{T}(undef, length(sol), 7)
    _fillputmatrix!(mat, sol, ind)
    return mat
end

function _fillputmatrix!(mat, sol, ind::Integer)
    for i ∈ eachindex(sol) 
        mat[i, 1] = sol[i][ind, 1]
        mat[i, 2] = sum(@view sol[i][ind, 2:5])
        mat[i, 3] = sum(@view sol[i][ind, 6:15])
        mat[i, 4] = sum(@view sol[i][ind, 16:21])
        mat[i, 5] = sum(@view sol[i][ind, 22:25])
        mat[i, 6] = sum(@view sol[i][ind, 26:27])
        mat[i, 7] = sol[i][ind, 28]
    end
end

function _fillputmatrix!(mat, sol, inds)
    for i ∈ eachindex(sol) 
        mat[i, 1] = sum(@view sol[i][inds, 1])
        mat[i, 2] = sum(@view sol[i][inds, 2:5])
        mat[i, 3] = sum(@view sol[i][inds, 6:15])
        mat[i, 4] = sum(@view sol[i][inds, 16:21])
        mat[i, 5] = sum(@view sol[i][inds, 22:25])
        mat[i, 6] = sum(@view sol[i][inds, 26:27])
        mat[i, 7] = sum(@view sol[i][inds, 28])
    end
end

makeoutputmatrix_susceptible(sol) = _makeoutputmatrix(sol, 1)
makeoutputmatrix_infectious(sol) = _makeoutputmatrix(sol, 2)
makeoutputmatrix_resistant(sol) = _makeoutputmatrix(sol, 3:5)
makeoutputmatrix_N(sol) = _makeoutputmatrix(sol, 1:5)
makeoutputmatrix_cumcases(sol) = _makeoutputmatrix(sol, 6)

function makeoutputmatrix_incidentcases(sol)
    incidentcases = makeoutputmatrix_cumcases(sol)
    ℓ = size(incidentcases, 1)
    for i ∈ axes(incidentcases, 1), j ∈ axes(incidentcases, 2)
        k = ℓ - i + 1
        if k == 1 
            incidentcases[k, j] = zero(incidentcases[k, j])
        else
            incidentcases[k, j] += -incidentcases[(k - 1), j]
        end
    end
    return incidentcases
end

const contacts_alllocations = makecontactmatrix(datadir("exp_raw", "MUestimates_all_locations_2.csv"))
const contacts_schools = makecontactmatrix(datadir("exp_raw", "MUestimates_school_2.csv"))
const contacts_notschools = contacts_alllocations .- contacts_schools
@assert minimum(contacts_notschools) >= 0
const contacts_home = makecontactmatrix(datadir("exp_raw", "MUestimates_home_2.csv"))
const contacts_other = makecontactmatrix(datadir("exp_raw", "MUestimates_other_locations_2.csv")) 
const contacts_lockdownschoolopen = (contacts_home .+ contacts_other .+ contacts_schools) .* 0.8
const contacts_lockdownschoolclosed = (contacts_home .+ contacts_other) .* 0.8

const births::Vector{Float64} = [
    55_690,
    57_781,
    60_041,
    59_046,
    58_791,
    58_590,
    58_027,
    56_014,
    56_725,
    55_098,
    54_488,
    52_861,
    51_308,
    49_863,
    46_809,
    47_786,
    46_959,
]

const save_positions = ( false, false )
# callback to age children by one year 
function ageswitch!(integrator)
    for g ∈ 1:15, j ∈ 1:5 
        i = 16 - g 
        integrator.u[j, (i + 1)] += integrator.u[j, i]
        integrator.u[j, i] = 0.0 
    end
end

startholidays!(integrator) = integrator.p.c = contacts_notschools
endholidays!(integrator) = integrator.p.c = contacts_alllocations
lockdownschoolopen!(integrator) = integrator.p.c = contacts_lockdownschoolopen
lockdownschoolclosed!(integrator) = integrator.p.c = contacts_lockdownschoolclosed

function updatebirthrate!(integrator)
    i = round(Int, integrator.t - 2005, RoundDown) 
    integrator.p.ν = births[i] 
end

const agetimes = [ x + 0.581 for x ∈ 2007:2023 ]
const agecb = PresetTimeCallback(agetimes, ageswitch!; save_positions)

const holidayfractions = [ 0.104, 0.255, 0.389, 0.498, 0.759, 0.893, 0.973 ]
const holidaytimes = [ vec([ x + f for x ∈ 2006:2019, f ∈ holidayfractions ]); [ 2020.104, 2021.603, 2021.759, 2021.893, 2021.973 ]; vec([ x + f for x ∈ 2022:2023, f ∈ holidayfractions ]) ]
const holidaycb = PresetTimeCallback(holidaytimes, startholidays!; save_positions)

const schoolfractions = [ 0.016, 0.123, 0.304, 0.403, 0.619, 0.805, 0.904 ]
const schooltimes = [ vec([ x + f for x ∈ 2006:2019, f ∈ schoolfractions ]); [ 2020.016, 2021.619, 2021.805, 2021.904 ]; vec([ x + f for x ∈ 2022:2023, f ∈ schoolfractions ]) ]
const schoolcb = PresetTimeCallback(schooltimes, endholidays!; save_positions)

const lockdownschoolopentimes = [ 2020.219, 2020.628, 2020.805, 2020.904, 2021.203, 2021.304, 2021.403 ]
const lockdownschoolopencb = PresetTimeCallback(lockdownschoolopentimes, lockdownschoolopen!; save_positions)

const lockdownschoolclosedtimes = [ 2020.221, 2020.759, 2020.893, 2020.973, 2021.255, 2021.389, 2021.498 ]
const lockdownschoolclosedcb = PresetTimeCallback(lockdownschoolclosedtimes, lockdownschoolclosed!; save_positions)

const yearstime = [ x + 0.581 for x ∈ 2007:2022 ]
const yearscb = PresetTimeCallback(yearstime, updatebirthrate!; save_positions)

# for convenience, set cumulative cases to 0 at the start of the analysis 
function resetcumcases!(integrator)
    for i ∈ 1:28
        integrator.u[6, i] = 0.0 
    end
end
const resetcumcasescb = PresetTimeCallback(2006.581, resetcumcases!; save_positions)



cbs = CallbackSet(agecb, holidaycb, schoolcb, lockdownschoolopencb, lockdownschoolclosedcb, yearscb)

#    redcb = PresetTimeCallback(reductiontime, reducetransmission!; save_positions)
#    rescb = PresetTimeCallback(reductiontime + 1, restoretransmission!; save_positions)
#    cbs = CallbackSet(redcb, rescb)



## test run 

# simplified to give a constant population 

# mortality estimated from mid-point of age ranges for women in 2020 from 
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/mortalityratesqxprincipalprojectionscotland
mortality = [
    506.44,  # 0
    29.53,  # 1
    21.56,  # 2
    7.04,  # 3
    14.07,  # 4
    10.38,  # 5
    13.39,  # 6
    3.25,  # 7
    3.15,  # 8
    3.15,  # 9
    19.32,  # 10
    6.46,  # 11
    12.64,  # 12
    18.90,  # 13
    9.40,  # 14
    26.27,  # 17
    32.97,  # 22
    47.73,  # 27
    71.67,  # 32
    90.65,  # 37
    196.06,  # 42
    276.98,  # 47
    428.51,  # 52
    723.39,  # 57
    1042.47,  # 62
    1625.24,  # 67
    2761.07, # 72
    100_000 / 14  # life expetancy at 75 in UK is 14 years 
    # https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/articles/lifeexpectancycalculator/2019-06-07
] ./ 100_000

#p = SirrrsParameters(0.1, contacts_alllocations, 48.7, births[1], mortality, 2.0, 0.75)
#p = SirrrsParameters(0.347, contacts_alllocations, 48.7, births[1], mortality, 1.0, 0.913)
p = SirrrsParameters(0.2, contacts_alllocations, 48.7, births[1], mortality, 1.0, 0.913)


pops = [ 
    # https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
    0,
    54_700,
    54_100,
    52_600,
    51_800,
    52_800,
    53_700,
    56_300,
    57_400,
    59_100,
    58_700,
    59_500,
    61_500,
    63_000,
    65_400,
    328_700,
    339_100,
    310_300,
    317_200,
    384_800,
    405_100,
    377_900,
    335_400,
    344_800,
    279_600,
    243_400,
    212_400,
    168_500 + 118_600 + 62_800 + 32_400
]

u0 = let 
    props = [ [ 0.3 ]; ones(20) .* 0.01; ones(4) .* 0.005; ones(2) .* 0.007; [ 0.0125 ] ]
    u = zeros(6, 28)
    for a ∈ 1:28
        u[1, a] = props[a] * pops[a] 
        u[2, a] = 0.01 * pops[a] 
        for i ∈ 3:5 
            u[i, a] = (pops[a] - props[a] * pops[a] - 0.01 * pops[a]) / 3
        end
    end
    u
end

tspan = ( 2006.581, 2023.5233 )
const saveat = collect(2016.7568:0.019165:2023.5233)

prob = ODEProblem(sirrrs!, u0, tspan, p)

# test with Memoization 

using BenchmarkTools 

solver_v2(prob, solver; kwargs...) = @memoize solve(prob, solver; kwargs...)


# test run
sol = solve(
    prob, Vern9(lazy=false); 
    p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=5e4, 
)

#=
julia> @benchmark sol = solve(prob, Vern9(lazy=false);
           p, u0, callback=cbs, saveat,
           abstol=1e-15, maxiters=5e4,
       )
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.059 s …   1.119 s  ┊ GC (min … max): 14.73% … 18.90%
 Time  (median):     1.090 s              ┊ GC (median):    14.30%
 Time  (mean ± σ):   1.090 s ± 22.901 ms  ┊ GC (mean ± σ):  15.83% ±  2.01%

  █                 █          █           █              █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.06 s         Histogram: frequency by time        1.12 s <

 Memory estimate: 1.27 GiB, allocs estimate: 13341552.
=#

#=
julia> @benchmark sol = solve(prob, Vern9(lazy=false);
                  p, u0, callback=cbs, saveat,
                  abstol=1e-15, maxiters=5e4,
              )
BenchmarkTools.Trial: 8 samples with 1 evaluation.
 Range (min … max):  541.239 ms … 691.330 ms  ┊ GC (min … max): 11.07% … 19.27%
 Time  (median):     611.782 ms               ┊ GC (median):    19.84%
 Time  (mean ± σ):   614.897 ms ±  50.827 ms  ┊ GC (mean ± σ):  16.95% ±  4.58%

  ▁              ▁▁       ▁       █                       ▁   ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁██▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁█ ▁
  541 ms           Histogram: frequency by time          691 ms <

 Memory estimate: 653.47 MiB, allocs estimate: 5547678.
=#

#=
@benchmark sol = solve(prob, Vern9(lazy=false);
                  p, u0, callback=cbs, saveat,
                  abstol=1e-15, maxiters=5e4,
              )
BenchmarkTools.Trial: 49 samples with 1 evaluation.
 Range (min … max):   82.406 ms … 151.593 ms  ┊ GC (min … max): 0.00% … 35.99%
 Time  (median):      99.555 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   103.017 ms ±  16.363 ms  ┊ GC (mean ± σ):  5.39% ± 11.15%

                ▂█   
  ▃▁▄▁▁▄▁▁▁▁▁▁▄▇██▄▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▃▁▄ ▁
  82.4 ms          Histogram: frequency by time          152 ms <

 Memory estimate: 46.00 MiB, allocs estimate: 294696.
=#

#=
julia> @benchmark solve(prob, Vern9(lazy=false);
           p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=5e4, save_idxs=[ 6 * x for x ∈ 1:28 ],
       )
BenchmarkTools.Trial: 48 samples with 1 evaluation.
 Range (min … max):   81.150 ms … 155.690 ms  ┊ GC (min … max): 0.00% … 39.40%
 Time  (median):      98.318 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   102.335 ms ±  16.122 ms  ┊ GC (mean ± σ):  5.99% ± 12.33%

              ▁▄█   
  ▃▁▁▁▃▁▁▁▃▄▁▇████▃▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▄▁▁▁▁▁▁▃ ▁
  81.2 ms          Histogram: frequency by time          156 ms <

 Memory estimate: 46.13 MiB, allocs estimate: 308819.
=#

#=
sol_out = makeoutputmatrices(sol)

fig = Figure(; size=( 1000, 1000 ))
ax1 = [ Axis(fig[i, 1]) for i ∈ 1:8 ]
ax2 = [ Axis(fig[i, 2]) for i ∈ 1:8 ]
ax3 = [ Axis(fig[i, 3]) for i ∈ 1:8 ]
ax4 = [ Axis(fig[i, 4]) for i ∈ 1:8 ]
for i ∈ 1:7 
    lines!(ax1[i], sol_out.t, sol_out.susceptible[:, i] ./ sol_out.N[:, i])
    lines!(ax2[i], sol_out.t, sol_out.infectious[:, i] ./ sol_out.N[:, i])
    lines!(ax3[i], sol_out.t, sol_out.resistant[:, i] ./ sol_out.N[:, i])
    lines!(ax4[i], sol_out.t, sol_out.N[:, i])
end
lines!(ax1[8], sol_out.t, [ sum(sol_out.susceptible[j, :]) #=./ sum(sol_out.N[j, :])=# for j ∈ axes(sol_out.susceptible, 1) ])
lines!(ax2[8], sol_out.t, [ sum(sol_out.infectious[j, :]) #=./ sum(sol_out.N[j, :])=# for j ∈ axes(sol_out.infectious, 1) ])
lines!(ax3[8], sol_out.t, [ sum(sol_out.resistant[j, :]) #=./ sum(sol_out.N[j, :])=# for j ∈ axes(sol_out.resistant, 1) ])
lines!(ax4[8], sol_out.t, [ sum(sol_out.N[j, :]) for j ∈ axes(sol_out.N, 1) ])

@unpack cumcases = sol_out






#formataxis!([ ax1; ax2; ax3; ax4 ]; hidex=true, hidey=true)

fig
=#
intagedatamatrix = let 
    agedatamatrix = zeros(Int(size(agedata, 1) / 7), 7)
    for i ∈ axes(agedatamatrix, 1), j ∈ axes(agedatamatrix, 2)
        agedatamatrix[i, j] = agedata[(7 * (i - 1) + j), (j + 11)]
    end
    # this is per 100_000 -- convert to numbers using the population numbers used above 
    agedatamatrix[:, 1] .*= (55_690 / 100_000)
    agedatamatrix[:, 2] .*= (sum(@view pops[2:5]) / 100_000)
    agedatamatrix[:, 3] .*= (sum(@view pops[6:15]) / 100_000)
    agedatamatrix[:, 4] .*= (sum(@view pops[16:21]) / 100_000)
    agedatamatrix[:, 5] .*= (sum(@view pops[22:25]) / 100_000)
    agedatamatrix[:, 6] .*= (sum(@view pops[26:27]) / 100_000)
    agedatamatrix[:, 7] .*= (pops[28] / 100_000)
    for k ∈ axes(agedatamatrix, 1) 
        i = size(agedatamatrix, 1) + 1 - k 
        i == 1 && continue
        if agedatamatrix[i, 1] > agedatamatrix[(i - 1), 1]
            for j ∈ axes(agedatamatrix, 2)
                agedatamatrix[i, j] += -agedatamatrix[(i - 1), j]
            end
        end
    end
    @assert minimum(agedatamatrix) >= 0

    round.(Int, agedatamatrix)
end




@model function fitmodel(data, prob, u0, cbs, mortality)
    #τ ~ truncated(Exponential(0.1), 1e-6, 10.0)
    τ ~ truncated(Exponential(0.1), 0.0, 10.0)
    ψ ~ truncated(Exponential(1), 0.0, 100.0) 
    ϕ ~ Beta(1, 1)

    p = SirrrsParameters(τ, contacts_alllocations, 48.7, births[1], mortality, ψ, 0.913)
    #@memoize solve(prob, Vern9(lazy=false); 
    #    p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=5e4,
    #)
    sol = solver_v2(
        prob, Vern9(lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=5e4,
    )
    #println("sol.retcode=$(sol.retcode), $(sol.retcode != :Success) when τ=$τ, ψ=$ψ, ϕ=$ϕ")
    if sol.retcode != :Success
        @info "Adding logprob -Inf when τ=$τ, ψ=$ψ, ϕ=$ϕ, sol.retcode=$(sol.retcode)"
        Turing.@addlogprob! -Inf
        return nothing
    end
    
    incidentcases = makeoutputmatrix_incidentcases(sol)

    for i ∈ axes(incidentcases, 1), j ∈ axes(incidentcases, 2)
        #data[i, j] ~ Binomial(round(Int, cumcases[i, j]), ϕ)
        #println("$(data[i, j]) ~ Binomial($(round(Int, cumcases[i, j])), $ϕ) $(data[i, j] ~ Binomial(round(Int, cumcases[i, j]), ϕ))")
        #println("logpdf(Binomial($(round(Int, cumcases[i, j])), $ϕ), $(data[i, j])) = $(logpdf(Binomial(round(Int, cumcases[i, j]), ϕ), data[i, j]))")
        #println("$(data[i, j]) ~ Poisson($(cumcases[i, j] * ϕ)) $(data[i, j] ~ Poisson(cumcases[i, j] * ϕ))")
        #println("logpdf(Poisson($(cumcases[i, j] * ϕ + 1e-10)), $(data[i, j])) = $(logpdf(Poisson(cumcases[i, j] * ϕ + 1e-10), data[i, j]))")
        data[i, j] ~ Poisson(incidentcases[i, j] * ϕ + 1e-10)
    end
end

#fitmodelp = fitmodel(agedatamatrix, u0, cbs, contacts_alllocations, births, mortality)

function fitmodel_target(data=intagedatamatrix, prob=prob, u0=u0, cbs=cbs, mortality=mortality)
    return Pigeons.TuringLogPotential(fitmodel(data, prob, u0, cbs, mortality))
end

const FitmodelType = typeof(fitmodel_target())
#=
p = SirrrsParameters(0.2, contacts_alllocations, 48.7, births[1], mortality, 1.0, 0.913)
sol = solve(prob, Vern9(lazy = false); 
    p, u0, callback = cbs, saveat, 
    abstol = 1e-15, reltol = 1e-15, maxiters = 5e6
)
=#
function Pigeons.initialization(target::FitmodelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext())
    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    Pigeons.update_state!(result, :τ, 1, 0.2)
    Pigeons.update_state!(result, :ψ, 1, 1.0)
    Pigeons.update_state!(result, :ϕ, 1, 0.005)

    return result
end
#=
@unpack cumcases = makeoutputmatrix(sol)
for i ∈ axes(cumcases, 1), j ∈ axes(cumcases, 2)
    ℓ = size(cumcases, 1)
    k = ℓ - i + 1
    if k == 1 
        cumcases[k, j] = 0.0 
    else
        cumcases[k, j] += -cumcases[(k - 1), j]
    end
end

sampledata = [
    rand(Binomial(round(Int, cumcases[i, j]), 0.01))
    for i ∈ axes(cumcases, 1), j ∈ axes(cumcases, 2)
]
=#
fitted_pt = pigeons( ;
    target=fitmodel_target(), 
    n_rounds=0,
    n_chains=16,
    multithreaded=true,
    record=[ traces; record_default() ],
    #seed=(1),
    variational=GaussianReference(),
)

@assert Pigeons.variable(fitted_pt.replicas[1].state, :τ) == [ 0.2 ]
@assert Pigeons.variable(fitted_pt.replicas[1].state, :ψ) == [ 1 ]


pt = increment_n_rounds!(fitted_pt, 1)
new_pt = pigeons(pt)
new_chains = Chains(new_pt)

resultdict = Dict(
    "chain" => new_chains, 
    "pt" => new_pt, 
    "n_rounds" => 1, 
    "n_chains" => 16,
)

safesave(datadir("sims", "firstchain.jld2"), resultdict)

pt = increment_n_rounds!(new_pt, 1)
new_pt = pigeons(pt)
new_chains = Chains(new_pt)

resultdict = Dict(
    "chain" => new_chains, 
    "pt" => new_pt, 
    "n_rounds" => 2, 
    "n_chains" => 16,
)

safesave(datadir("sims", "firstchain.jld2"), resultdict)

pt = increment_n_rounds!(new_pt, 1)
new_pt = pigeons(pt)
new_chains = Chains(new_pt)

resultdict = Dict(
    "chain" => new_chains, 
    "pt" => new_pt, 
    "n_rounds" => 3, 
    "n_chains" => 16,
)

safesave(datadir("sims", "firstchain.jld2"), resultdict)

pt = increment_n_rounds!(new_pt, 1)
new_pt = pigeons(pt)
new_chains = Chains(new_pt)

resultdict = Dict(
    "chain" => new_chains, 
    "pt" => new_pt, 
    "n_rounds" => 4, 
    "n_chains" => 16,
)

safesave(datadir("sims", "firstchain.jld2"), resultdict)

###
###
p = SirrrsParameters(7.0697, contacts_alllocations, 48.7, births[1], mortality, 5.1813, 0.913)

prob = ODEProblem(sirrrs!, u0, tspan, p)

# test run
sol = solve(prob, Vern9(lazy=false); 
    p, u0, callback=cbs, saveat, 
    abstol=1e-15, maxiters=5e4,
)


@unpack t, incidentcases = makeoutputmatrices(sol)

fig = Figure(; size=( 1000, 1000 ))
ax1 = [ Axis(fig[i, 1]) for i ∈ 1:8 ]
for i ∈ 1:7 
    lines!(ax1[i], t, incidentcases[:, i] .* 0.0545)
    scatter!(ax1[i], t, intagedatamatrix[:, i]; color=:black, markersize=3)
end
lines!(ax2[8], t, [ sum(incidentcases[j, :]) .* 0.0545 #=./ sum(sol_out.N[j, :])=# for j ∈ axes(sol_out.infectious, 1) ])


#formataxis!([ ax1; ax2; ax3; ax4 ]; hidex=true, hidey=true)

fig
=#

#=
using ModelingToolkit
de = modelingtoolkitize(prob)
ModelingToolkit.generate_jacobian(de)[2] 
prob_jac2 = ODEProblem(de, [], tspan; jac=true)

sol = solve(prob_jac2, Vern9(lazy=false); 
    p, u0, callback=cbs, saveat, 
    abstol=1e-15, maxiters=5e4,
)
=#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple models without having a large 
# number of values in global scope potentially leading to unexpected results
rsvparms = let 
    γ = 48.7    # generation time 7.5 days
    μ = .0087   # Scotland's birth rate = 48000 / 5.5e6
    ϕ = -.5π

    # When is the infection parameter expected to change?
    # Find dates (as fractions of year) when Strigency Index goes above then below 50
    inds = findall(x -> x >= 50, crgtdata.StringencyIndex_Average)
    reduceday = crgtdata.Date[inds[1]]
    increaseday = crgtdata.Date[last(inds)]
    println("Stringency ≥ 50 on $(printrawdate(crgtdata.RawDate[inds[1]]))")
    println("Stringency < 50 on $(printrawdate(crgtdata.RawDate[last(inds)]))")
    # note the Stringency Index is plotted by code in `npisimulation.jl`

    # Vector of recorded cases 
    casesvector = data.Cases 

    ## Times when we have data pre-lockdown (i.e. times to save simulation)
    # add one pre-data date so that we can calculate a weekly incidence for the first data point
    savetimes = [ minimum(data.Date) - 7 / 365; data.Date ] # 354 elements

    @ntuple ϕ γ μ casesvector reduceday increaseday savetimes
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Optimize magnitude of the effect of non-pharmaceutical interventions 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Callbacks for optimization 
opt_cbs = let 
    @unpack increaseday, reduceday = rsvparms
    save_positions = ( false, false )
    resetcb = PresetTimeCallback(reduceday + 1e-9, opt_restoretransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, opt_changetransmission!; save_positions)
    CallbackSet(resetcb, rescb)
end

multipliers_psi0 = let 
    @unpack casesvector, increaseday, reduceday, savetimes = rsvparms
    R0 = 1.215
    β1 = .1
    ψ = 0
    initialvalues = [ .8, 1. ]
    iterativeopt(R0, β1, ψ, initialvalues; 
        casesvector, increaseday, opt_cbs, reduceday, savetimes) 
end

multipliers_psi5 = let 
    @unpack casesvector, increaseday, reduceday, savetimes = rsvparms
    R0 = 1.285
    β1 = .082
    ψ = 5
    initialvalues = [ .8, 1. ]
    iterativeopt(R0, β1, ψ, initialvalues;
        casesvector, increaseday, opt_cbs, reduceday, savetimes) 
end

multipliers_psi13_2 = let 
    @unpack casesvector, increaseday, reduceday, savetimes = rsvparms
    R0 = 1.6
    β1 = 0
    ψ = 13.2
    initialvalues = [ .8, 1. ]
    iterativeopt(R0, β1, ψ, initialvalues; 
        casesvector, increaseday, opt_cbs, reduceday, savetimes) 
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define the ODE problem 
prob = let 
    # different parameters may be used in the simulations
    @unpack ϕ, γ, μ = rsvparms
    R0 = 1.6 
    β0 = R0 * (γ + μ)
    ω = 365.25 / 400
    β1 = .0 
    ψ = 0

    tspan = ( 1015.35, 2023.6 )
    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, 1., 1.) 
    I0 = .007
    S0 = .5
    u0 = sirns_u0(S0, I0; equalrs = true, p, t0 = 1015.35)
    ODEProblem(sirns!, u0, tspan, p)
end

## Run the simulations 

# Without immune boosting
rsvsim_psi0 = let 
    @unpack ϕ, γ, μ, reduceday, increaseday, savetimes = rsvparms
    βreduction, βreturn = multipliers_psi0 
    R0 = 1.215
    ψ = 0
    immuneduration = .5
    β1 = .1
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    
    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = .35)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)
    
    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 1e7)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

rsvsim_psi5 = let 
    @unpack ϕ, γ, μ, reduceday, increaseday, savetimes = rsvparms
    βreduction, βreturn = multipliers_psi5
    R0 = 1.285
    ψ = 5
    immuneduration = .5
    β1 = .082
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = .35) 
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 1e7)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 

rsvsim_psi13_2 = let 
    @unpack ϕ, γ, μ, reduceday, increaseday, savetimes = rsvparms
    βreduction, βreturn = multipliers_psi13_2
    R0 = 1.6
    ψ = 13.2
    immuneduration = .5
    β1 = ϕ = 0
    θ = .0025
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, βreduction, βreturn) 
    u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = .35)
    save_positions = ( false, false )
    redcb = PresetTimeCallback(reduceday, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    cbs = CallbackSet(redcb, rescb)

    sol = solve(prob, Vern9(lazy = false); 
        p, u0, callback = cbs, saveat = savetimes, 
        abstol = 1e-15, reltol = 1e-15, maxiters = 5e6)
    compartments = modelcompartments(sol, p)
    cases = casespertimeblock(compartments[:cc]) * 5_500_000 * θ
    @ntuple compartments cases
end 
