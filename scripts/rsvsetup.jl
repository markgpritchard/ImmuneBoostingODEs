
# This script is called by `rsvanalysis.jl` and `displayrsvanalysis.jl`


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

## Times to save simulations
saveat = let 
    savefirst = data.Date[1] - 7 / 365  # to allow calculation of new cases in first week
    [ [ savefirst ]; data.Date ]
end

## Callbacks  
cbs = let 
    # When is the infection parameter expected to change?
    # Find dates (as fractions of year) when Strigency Index goes above then below 50
    inds = findall(x -> x >= 50, crgtdata.StringencyIndex_Average)
    reduceday = crgtdata.Date[inds[1]]
    increaseday = crgtdata.Date[last(inds)]
    save_positions = ( false, false )
    resetcb = PresetTimeCallback(reduceday + 1e-9, reducetransmission!; save_positions)
    rescb = PresetTimeCallback(increaseday, restoretransmission!; save_positions)
    CallbackSet(resetcb, rescb)
end
