# This script is called by `rsvanalysis.jl` and `displayrsvanalysis.jl`

# RSV data from Scotland
data = processrsvdata("respiratory_scot.csv", "rsv.csv")

# Age-specific data
agedata = processagedata("respiratory_age.csv", "rsv_age.csv")

# Mobility data
mobilitydata = processmobilitydata(
    "2020_GB_Region_Mobility_Report.csv", 
    "2021_GB_Region_Mobility_Report.csv", 
    "2022_GB_Region_Mobility_Report.csv"
)

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

## Times to save simulations
saveat = let 
    savefirst = data.Date[1] - 7 / 365  # to allow calculation of new cases in first week
    [ [ savefirst ]; data.Date ]
end

## Callbacks  
const MOBILITYCALLBACKTIMES = mobilitydata.gtdate
const MOBILITYCALLBACKPROPORTIONS = mobilitydata.reduction

function mobilityaffect!(integrator)
    index = findfirst(x -> x >= integrator.t, MOBILITYCALLBACKTIMES)
    reduction = 1 - MOBILITYCALLBACKPROPORTIONS[index] 
    adjustedreduction = integrator.p.betaprimemultiplier * reduction
    newbetazero = *(
        integrator.p.originalβ0,
        (1 - adjustedreduction)
    ) 
    newparameters = SirnsParameters(
        newbetazero, 
        integrator.p.β1, 
        integrator.p.ϕ, 
        integrator.p.γ, 
        integrator.p.μ, 
        integrator.p.ψ, 
        integrator.p.ω, 
        integrator.p.originalβ0, 
        integrator.p.betaprimemultiplier,
        integrator.p.finalbetaprime, 
        integrator.p.proportiondetected
    )
    integrator.p = newparameters
end

function finalmobilityaffect!(integrator)
    newbetazero = *(
        integrator.p.originalβ0, 
        integrator.p.finalbetaprime, 
    ) 
    newparameters = SirnsParameters(
        newbetazero, 
        integrator.p.β1, 
        integrator.p.ϕ, 
        integrator.p.γ, 
        integrator.p.μ, 
        integrator.p.ψ, 
        integrator.p.ω, 
        integrator.p.originalβ0, 
        integrator.p.betaprimemultiplier,
        integrator.p.finalbetaprime, 
        integrator.p.proportiondetected
    )
    integrator.p = newparameters
end

cbs = let 
    save_positions = ( false, false )
    mobilitycb = PresetTimeCallback(MOBILITYCALLBACKTIMES, mobilityaffect!; save_positions)
    finalmobilitycb = PresetTimeCallback(
        last(MOBILITYCALLBACKTIMES) + 1 / 365, finalmobilityaffect!; 
        save_positions
    )
    CallbackSet(mobilitycb, finalmobilitycb) 
end

optimearlyaffect!(integrator) = integrator.p[7] = 1.0

function optimmobilityaffect!(integrator)
    index = findfirst(x -> x >= integrator.t, MOBILITYCALLBACKTIMES)
    reduction = 1 - MOBILITYCALLBACKPROPORTIONS[index] 
    adjustedreduction = ImmuneBoostingODEs._logistic(integrator.p[5]) * reduction
    integrator.p[8] = max(1 - adjustedreduction, 0.0)  # never negative
end

function optimfinalmobilityaffect!(integrator)
    integrator.p[8] = ImmuneBoostingODEs._logistic(integrator.p[6])
end

optimcbs = let 
    save_positions = ( false, false )
    # reset `betaprime` which gets mutated in each iteration
    optimearlycallback = PresetTimeCallback(
        1996.737, optimearlyaffect!; 
        save_positions
    )
    mobilitycb = PresetTimeCallback(
        MOBILITYCALLBACKTIMES, optimmobilityaffect!; 
        save_positions
    )
    finalmobilitycb = PresetTimeCallback(
        last(MOBILITYCALLBACKTIMES) + 1 / 365, optimfinalmobilityaffect!; 
        save_positions
    )
    CallbackSet(optimearlycallback, mobilitycb, finalmobilitycb) 
end

