
using DataFrames, DifferentialEquations, Random, Turing

include("samplepriors.jl")

poissoncases(::Missing) = missing 
poissoncases(x::Number) = rand(Poisson(x))

Random.seed!(1729)

if isfile(datadir("sims", "priormodeloutputs.jld2"))
    @info "Loading prior model outputs"
    priormodeloutputs = load(datadir("sims", "priormodeloutputs.jld2"))
else
    priormodeloutputs = let 
        pmo = Dict{String, Matrix{Union{Missing, Float64}}}()
        names = [ "priors$x" for x ∈ [ "01", "02", "04", "1", "2", "4", "6" ] ]
        omegas = [ 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 6.0 ]
        for i ∈ 1:7 
            priorsdf = DataFrame(priorsdict[names[i]])
            push!(pmo, names[i] => runfittedsimulations(priorsdf, omegas[i], saveat, cbs))
        end
        pmo
    end
    safesave(datadir("sims", "priormodeloutputs.jld2"), priormodeloutputs)
end

if isfile(datadir("sims", "priormodelcases.jld2"))
    @info "Loading prior model number of cases"
    priormodelcases = load(datadir("sims", "priormodelcases.jld2"))
else
    priormodelcases = let 
        pmc = Dict{String, Matrix{Union{Missing, Int}}}()
        names = [ "priors$x" for x ∈ [ "01", "02", "04", "1", "2", "4", "6" ] ]
        omegas = [ 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 6.0 ]
        Threads.@threads for i ∈ 1:7 
            push!(pmc, names[i] => poissoncases.(priormodeloutputs[names[i]]))
        end
        pmc
    end
    safesave(datadir("sims", "priormodelcases.jld2"), priormodelcases)
end

priorsvaluequantiles = Vector{Vector{Vector{Float64}}}(undef, 7)
Threads.@threads for i ∈ 1:7
    priorsvaluequantiles[i] = fittedsimulationquantiles(
        priormodeloutputs["priors$([ "01", "02", "04", "1", "2", "4", "6" ][i])"],
        [ 0.025, 0.5, 0.975 ]
    )
end

priorscasesquantiles = Vector{Vector{Vector{Float64}}}(undef, 7)
Threads.@threads for i ∈ 1:7
    priorscasesquantiles[i] = fittedsimulationquantiles(
        priormodelcases["priors$([ "01", "02", "04", "1", "2", "4", "6" ][i])"], 
        [ 0.025, 0.5, 0.975 ]
    )
end
