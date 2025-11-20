
# This script is called by `rsvanalysis.jl` and `displayrsvanalysis.jl`

# RSV data from Scotland
data = processrsvdata("rsv.csv", "respiratory_scot.csv")

# Age-specific data
agedata = processagedata("rsv_age.csv", "respiratory_age.csv")

# Data from Oxford Covid-19 Government Response Tracker
crgtdata = processcrgtvdata("crgt.csv", "OxCGRT_compact_subnational_v1.csv")

# Google mobility data (as a const as it is called by ODE callbacks)
const MOBILITYDATA = processmobilitydata(
    "mobility.csv",
    "2020_GB_Region_Mobility_Report.csv",
    "2021_GB_Region_Mobility_Report.csv",
    "2022_GB_Region_Mobility_Report.csv",
)

# To avoid splitting outbreaks, count cases from April each year 
let 
    april1value = 90 / 365  # NB assuming 28 days in February
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
    [[savefirst]; data.Date]
end

## Callbacks  

function updatereduction!(integrator)
    ind = findfirst(x -> x >= integrator.t, MOBILITYDATA.fractiondate)
    integrator.p = SirnsParameters(integrator.p, MOBILITYDATA.betareduction[ind])
    return nothing
end

betareductioncallback = PresetTimeCallback(
    MOBILITYDATA.fractiondate, updatereduction!; 
    save_positions=(false, false),  # `fitmodel` assumes saving only at specified `saveat` times
)
