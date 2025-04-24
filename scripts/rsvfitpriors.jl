
include("rsvfitmodel.jl")

Random.seed!(1729)

if isfile(datadir("sims", "priorsdict.jdl2"))
    @info "Loading prior values"
    priorsdict = load(datadir("sims", "priorsdict.jdl2"))
else
    priorsdict = let 
        prob = fittedsimulationsetup(saveat)
        model01 = fitmodel(data.Cases, prob, cbs, saveat; omega=0.1)
        priors01 = sample(model01, Prior(), MCMCThreads(), 4000, 4)
        model02 = fitmodel(data.Cases, prob, cbs, saveat; omega=0.2)
        priors02 = sample(model02, Prior(), MCMCThreads(), 4000, 4)
        model04 = fitmodel(data.Cases, prob, cbs, saveat; omega=0.4)
        priors04 = sample(model04, Prior(), MCMCThreads(), 4000, 4)
        model1 = fitmodel(data.Cases, prob, cbs, saveat; omega=1.0)
        priors1 = sample(model1, Prior(), MCMCThreads(), 4000, 4)
        model2 = fitmodel(data.Cases, prob, cbs, saveat; omega=2.0)
        priors2 = sample(model2, Prior(), MCMCThreads(), 4000, 4)
        model4 = fitmodel(data.Cases, prob, cbs, saveat; omega=4.0)
        priors4 = sample(model4, Prior(), MCMCThreads(), 4000, 4)
        model6 = fitmodel(data.Cases, prob, cbs, saveat; omega=6.0)
        priors6 = sample(model6, Prior(), MCMCThreads(), 4000, 4)

        @strdict priors01 priors02 priors04 priors1 priors2 priors4 priors6
    end

    safesave(datadir("sims", "priorsdict.jld2"), priorsdict)
end

priorsvector = [
    priorsdict["priors01"],
    priorsdict["priors02"],
    priorsdict["priors04"],
    priorsdict["priors1"],
    priorsdict["priors2"],
    priorsdict["priors4"],
    priorsdict["priors6"],    
]

priorsvalues = [
    fittedsimulationquantiles(DataFrame(p), omega, saveat, cbs)
    for (p, omega) ∈ zip(priorsvector, [ 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 6.0 ])
]
