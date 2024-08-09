
using DrWatson
@quickactivate :ImmuneBoostingODEs

using CairoMakie, DataFrames, DifferentialEquations, Pigeons, Random, Turing


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RSV data from Scotland
data = processrsvdata("respiratory_scot.csv", "rsv.csv")

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

saveat = let  # copy of list in `rsvanalysis.jl`
    savefirst = data.Date[1] - 7 / 365  # to allow calculation of new cases in first week
    [ [ savefirst ]; data.Date ]
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load results 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maxrounds = 12 

for i ∈ 1:4 
    _loaded = false 
    j = maxrounds
    while !_loaded && j >= 1
        filename = "rsvparameters_omega_0.1_id_$(i)_nrounds_$(j).jld2"
        if isfile(datadir("sims", filename))
            if i == 1 
                global rsvparameters01 = DataFrame(load(datadir("sims", filename))["chain"])
            else 
                _rsvparameters01 = DataFrame(load(datadir("sims", filename))["chain"])
                _rsvparameters01.chain = [ i for _ ∈ axes(_rsvparameters01, 1) ]
                global rsvparameters01 = vcat(rsvparameters01, _rsvparameters01)
            end 
            _loaded = true
        end
        j += -1
    end
end

plotchains(rsvparameters01)

function loadrsvdata(omega; maxrounds=12)
    df = DataFrame(
        :iteration => Int[ ],
        :chain => Int[ ],
        :β0 => Float64[ ],
        :β1 => Float64[ ],
        :ϕ => Float64[ ],
        :ψ => Float64[ ],
        :βreduction1 => Float64[ ],
        :βreduction2 => Float64[ ],
        :detection => Float64[ ],
        :log_density => Float64[ ],
    )
    for i ∈ 1:4 
        _loaded = false 
        j = maxrounds
        while !_loaded && j >= 1
            filename = "rsvparameters_omega_$(omega)_id_$(i)_nrounds_$(j).jld2"
            if isfile(datadir("sims", filename))
                _df = DataFrame(load(datadir("sims", filename))["chain"])
                _df.chain = [ i for _ ∈ axes(_df, 1) ]
                df = vcat(df, _df)
                _loaded = true
            end
            j += -1
        end
    end
    return df
end


rsvparameters01 = loadrsvdata(0.1)
plotchains(rsvparameters01)

rsvparameters025 = loadrsvdata(0.25)
plotchains(rsvparameters025)

rsvparameters1 = loadrsvdata(1.0)
plotchains(rsvparameters1)

rsvparameters2 = loadrsvdata(2.0)
plotchains(rsvparameters2)

rsvparameters4 = loadrsvdata(4.0)
plotchains(rsvparameters4)

rsvparameters26 = loadrsvdata(26.0)
plotchains(rsvparameters26)
