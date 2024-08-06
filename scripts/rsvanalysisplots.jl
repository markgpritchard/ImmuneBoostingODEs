
using DrWatson
@quickactivate :ImmuneBoostingODEs
using CairoMakie, DataFrames, DifferentialEquations, Pigeons, Turing

maxruns = 10

for id ∈ 1:4 
    for i ∈ 1:maxruns 
        nextfilename = "rsvparameters_id_$(id)_nrounds_$(i + 1).jld2"
        isfile(datadir("sims", nextfilename)) && continue
        filename = "rsvparameters_id_$(id)_nrounds_$(i).jld2"
        if isfile(datadir("sims", filename))
            if id == 1
                global samplevalues = DataFrame(load(datadir("sims", filename))["chain"])
            else
                sv = DataFrame(load(datadir("sims", filename))["chain"])
                sv.chain = [ id for _ ∈ axes(sv, 1) ]
                append!(samplevalues, sv) 
            end
        else 
            continue
        end
    end
end

include("rsvanalysisconsts.jl")


plotchains(samplevalues)

p = SirrrsParameters(0.2, contacts_alllocations, 48.7, births[1], mortality, 1.0, 0.913)

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

prob = ODEProblem(sirrrs!, u0, tspan, p)

sampleincidence = solvesamples(
    prob, samplevalues; 
    births, callbackset=cbs, contacts_alllocations, mortality, saveat, u0
)

fig = Figure(; size=( 800, 800 ))
axs = [ Axis(fig[i, 1]) for i ∈ 1:7 ]
for i ∈ 1:7 
    scatter!(axs[i], saveat[2:end], intagedatamatrix[2:end, i]; color=:black, markersize=3)
    lines!(axs[i], saveat[2:end], sampleincidence.means[2:end, i])
    band!(
        axs[i], saveat[2:end], sampleincidence.lq[2:end, i], sampleincidence.uq[2:end, i]; 
        color=( :gray, 0.25 )
    )

end

fig
