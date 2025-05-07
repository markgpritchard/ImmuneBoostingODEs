
using DrWatson

@quickactivate :ImmuneBoostingODEs

#using AdvancedHMC, DataFrames, DifferentialEquations, Optim, Pathfinder, Random, Turing
#using AdvancedHMC, DataFrames, DifferentialEquations, Optim, Random, Turing
#using AdvancedHMC, DataFrames, DifferentialEquations, Pathfinder, Random, Turing
using DataFrames, DifferentialEquations, Optimization, OptimizationOptimJL, Random, Turing 

#using Transducers: ThreadedEx

include("samplepriors.jl")

testrun = true 

if length(ARGS) == 2 
    omega = parse(Float64, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
else
    omega = 2.0
    n_rounds = testrun ? 25 : 10_000
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prob = fittedsimulationsetup(saveat)

include("rsvfitmodel.jl")

#const constcases = data.Cases
#=
transformtheta(x) = exp(x) / (1 + exp(x)) =#

function transformforoptim(v)
    @assert length(v) == 5 
    return [
        _transformforoptim(v[1] / 50),
        _transformforoptim(v[2]),
        _transformforoptim(v[3] / 2π + 0.5),
        v[4],
        _transformforoptim(v[5]),
    ]
end

_transformforoptim(x) = log(x / (1 - x))

function transformintooptim(v)
    @assert length(v) == 5 
    return [
        50 * _transformintooptim(v[1]),
        _transformintooptim(v[2]),
        2π * _transformintooptim(v[3]) - π,
        v[4],
        _transformintooptim(v[5]),
        ]
end

_transformintooptim(x) = exp(x) / (1 + exp(x))


function optimmodel(u, p)
    ω, data, prob, callback, saveat = p
    # take any real values in paramvector and transform into "acceptable" values 
    R0, β1, ϕ, ψt, reduce1  = transformintooptim(u)
    ψ = exp(0.7 * ψt)
    reduce2 = 0.5 + reduce1 / 2
    detection = 0.015

    β0 = R0 * 48.7087
    p = SirnsParameters(
        β0, 
        β1, 
        ϕ, 
        48.7,  # γ 
        0.0087,  # μ 
        ψ, 
        ω, 
        β0, 
        reduce1 * β0, 
        reduce2 * β0
    )
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, 
        u0, 
        callback, 
        saveat, 
        save_idxs=[ 8 ], 
        abstol=1e-15, 
        maxiters=1e8, 
        verbose=false,
    )
    #println(sol)
    if sol.retcode != :Success
        println("Fail: $(transformintooptim(u)) -> $Inf")
        return Inf
    end
    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000
    #incidentcases = round.(Int, max.(1, casespertimeblock(cumulativecases) .* 5_450_000))
    println("$(transformintooptim(u)) -> $(sum(cumulativecases)) cases -> $(sum([ abs2(data[i] - incidentcases[i] * detection) for i ∈ eachindex(data) ]))")
    return sum([ abs2(data[i] - incidentcases[i] * detection) for i ∈ eachindex(data) ])
    #println("$(transformintooptim(u)) -> $(sum(cumulativecases)) cases -> $(sum([ -log(pdf(NegativeBinomial(incidentcases[i] * 0.01 * 0.25 + 1e-10, 0.2), data[i])) for i in eachindex(data) ]))")
    #return sum([ -log(pdf(NegativeBinomial(incidentcases[i] * 0.01 * 0.25 + 1e-10, 0.2), data[i])) for i in eachindex(data) ])
#        #incidence[i] ~ Poisson(incidentcases[i] + 1e-10)

end


poissonerror = [ -log(pdf(Poisson(incidentcases[i] * 0.01 + 1e-10), data.Cases[i])) for i in eachindex(data.Cases) ]

μ = [ incidentcases[i] * 0.01 for i in eachindex(data.Cases) ]
r = m * 0.2 / 0.8

nberror = [ -log(pdf(NegativeBinomial(incidentcases[i] * 0.01 * 0.25, 0.2), data.Cases[i])) for i in eachindex(data.Cases) ]

v = m / p

#=

result1 = optimize(
    x -> optimmodel(x, omega, data.Cases; callback=cbs, saveat), 
    [ 1.0, -2.0, 0.0, 4.0, 1.0, 2.0, -4.5 ], 
    NelderMead(), 
    Optim.Options(; iterations=5000)
)
result2 = optimize(
    x -> optimmodel(x, omega, data.Cases; callback=cbs, saveat), 
    [ 2.0, -2.0, 0.0, -4.0, 1.0, 2.0, -4.5 ], 
    NelderMead(), 
    Optim.Options(; iterations=5000)
)
result3 = optimize(
    x -> optimmodel(x, omega, data.Cases; callback=cbs, saveat), 
    [ 3.0, -2.0, 0.0, 4.0, 1.0, 2.0, -4.5 ], 
    NelderMead(), 
    Optim.Options(; iterations=5000)
)
result4 = optimize(
    x -> optimmodel(x, omega, data.Cases; callback=cbs, saveat), 
    [ 10.0, -2.0, 0.0, -4.0, 1.0, 2.0, -4.5 ], 
    NelderMead(), 
    Optim.Options(; iterations=5000)
=#
#=
optimvalues = Vector{Vector{Float64}}(undef, 16) 
let 
    initvalues = [
        [ 1.0, 0.001, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.001, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 1.0, 0.001, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.001, 0.0, 3.0, 0.5, 0.75, 0.01 ],
        [ 1.0, 0.5, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.5, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 1.0, 0.5, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.5, 0.0, 3.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.001, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.001, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.001, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.001, 0.0, 3.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.5, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.5, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.5, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.5, 0.0, 3.0, 0.5, 0.75, 0.01 ],
    ]
    Threads.@threads for i ∈ 1:16
        r = optimize(
            x -> optimmodel(x, omega, data.Cases; callback=cbs, saveat), 
            transformforoptim(initvalues[i]), 
            NelderMead(), 
            Optim.Options(; iterations=optimiterations)
        )
        optimvalues[i] = transformintooptim(Optim.minimizer(r))
    end
end


pathfindervalues = multipathfinder(
    fitmodel(data.Cases, prob, cbs, saveat; omega), 
    1000; 
    executor=ThreadedEx(),
    init=[
        [ 1.0, 0.001, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.001, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 1.0, 0.001, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.001, 0.0, 3.0, 0.5, 0.75, 0.01 ],
        [ 1.0, 0.5, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.5, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 1.0, 0.5, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 1.0, 0.5, 0.0, 3.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.001, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.001, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.001, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.001, 0.0, 3.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.5, 0.0, -7.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.5, 0.0, -7.0, 0.5, 0.75, 0.01 ],
        [ 10.0, 0.5, 0.0, 3.0, 0.95, 0.975, 0.01 ],
        [ 10.0, 0.5, 0.0, 3.0, 0.5, 0.75, 0.01 ],
    ],
    maxtime=7200,
)
=#

#optimfunction = OptimizationFunction(optimmodel, AutoFiniteDiff())
optimfunction = OptimizationFunction(optimmodel, Optimization.AutoForwardDiff())
optimvalues = Vector{Vector{Float64}}(undef, 16)
let 
    initvalues = [
        [ 1.0, 0.001, 0.0, -7.0, 0.95 ],
        [ 1.0, 0.001, 0.0, -7.0, 0.5 ],
        [ 1.0, 0.001, 0.0, 3.0, 0.95 ],
        [ 1.0, 0.001, 0.0, 3.0, 0.5 ],
        [ 1.0, 0.5, 0.0, -7.0, 0.95 ],
        [ 1.0, 0.5, 0.0, -7.0, 0.5 ],
        [ 1.0, 0.5, 0.0, 3.0, 0.95 ],
        [ 1.0, 0.5, 0.0, 3.0, 0.5 ],
        [ 10.0, 0.001, 0.0, -7.0, 0.95 ],
        [ 10.0, 0.001, 0.0, -7.0, 0.5 ],
        [ 10.0, 0.001, 0.0, 3.0, 0.95 ],
        [ 10.0, 0.001, 0.0, 3.0, 0.5 ],
        [ 10.0, 0.5, 0.0, -7.0, 0.95 ],
        [ 10.0, 0.5, 0.0, -7.0, 0.5 ],
        [ 10.0, 0.5, 0.0, 3.0, 0.95 ],
        [ 10.0, 0.5, 0.0, 3.0, 0.5 ],
    ]
    Threads.@threads for i ∈ 1:16
        optimproblem = OptimizationProblem(
            optimfunction, 
            transformforoptim(initvalues[i]), 
            ( omega, data.Cases, prob, cbs, saveat )
        )
        optimsolution = solve(optimproblem, BFGS(); maxtime=900)
        optimvalues[i] = transformintooptim(optimsolution)
    end
end

optimproblem = OptimizationProblem(
            optimfunction, 
            transformforoptim([ 20.0, 0.001, 0.0, -7.0, 0.95 ]), 
            ( omega, data.Cases, prob, cbs, saveat )
        )
        optimsolution = solve(optimproblem, LBFGS(); maxtime=100)
        #optimsolution = solve(optimproblem, Optim.NelderMead(); maxtime=100)

#=
optimproblem = OptimizationProblem(
    optimfunction, 
    transformforoptim([ 1.0, 0.001, 0.0, -7.0, 0.95 ]), 
    ( omega, data.Cases, prob, cbs, saveat )
)
optimsolution = solve(optimproblem, BFGS())
=#
#=
initvalues1 = Optim.minimizer(result1)
initvalues2 = Optim.minimizer(result2)
initvalues3 = Optim.minimizer(result3)
initvalues4 = Optim.minimizer(result4)
=#
chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat; omega),
    Turing.NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    #16;
    #initial_params=optimvalues,
    initial_params=[
        [ 20.0, 0.164, -0.912, -7.0, 0.886, 0.933, 0.055, 0.0094 ],
        [ 1.0, 0.040, -0.124, -7.0, 0.783, 0.975, 0.086, 0.0075 ],
        [ 20.0, 0.032, 1.331, 3.0, 0.870, 0.883, 0.146, 0.0204 ],
        [ 1.0, 0.166, 0.378, 3.0, 0.629, 0.845, 0.089, 0.0138 ],
    ],
)

chaindf = DataFrame(chain)
plotchains(chaindf)


chain2 = sample(
    fitmodel(data.Cases, prob, cbs, saveat; omega),
    Turing.NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
)
chaindf2 = DataFrame(chain2)
plotchains(chaindf2)
sort!

#omegapriors = priorsdict["priors2"]
omegapriors = priorsdict["$omega"]
filteredomegapriorsdf = DataFrame(omegapriors)
filter!(:rzero_t => x -> x > 0, filteredomegapriorsdf)
sort!(filteredomegapriorsdf, :lp; rev=true)

function initialparamsdf(df, i)
    return [ 
        df.rzero_t[i],
        df.betaone_t[i], 
        df.phi_t[i], 
        df.psi_t[i], 
        df.betareduction1_t[i], 
        df.betareduction2_t[i], 
        df.pparameter_t[i], 
        df.detection_t[i] 
    ]
end

chain3 = sample(
    fitmodel(data.Cases, prob, cbs, saveat; omega),
    Turing.NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    initial_params=[
        initialparamsdf(filteredomegapriorsdf, 1),
        initialparamsdf(filteredomegapriorsdf, 2),
        initialparamsdf(filteredomegapriorsdf, 3),
        initialparamsdf(filteredomegapriorsdf, 4),
    ],
)
chaindf3 = DataFrame(chain3)
plotchains(chaindf3)


#=


chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat; omega),
    NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    init_params=collect.(eachrow(result_multi.draws_transformed.value[1:n_chains, :, 1])),
)


Random.seed!(round(Int, omega * 100))

chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat; omega),
    NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    initial_params=[
        [
            initvalues1[1],  # R0
            transformtheta(initvalues1[2]),  # β1
            2π * transformtheta(initvalues1[3]) - π,  # ϕ
            initvalues1[4],  # ψt
            transformtheta(initvalues1[5]),  # βreduction1
            transformtheta(initvalues1[6]),  # βreduction2
            transformtheta(initvalues1[7])  # detection
        ],
        [
            initvalues2[1],  # R0
            transformtheta(initvalues2[2]),  # β1
            2π * transformtheta(initvalues2[3]) - π,  # ϕ
            initvalues2[4],  # ψt
            transformtheta(initvalues2[5]),  # βreduction1
            transformtheta(initvalues2[6]),  # βreduction2
            transformtheta(initvalues2[7])  # detection
        ],
        [
            initvalues3[1],  # R0
            transformtheta(initvalues3[2]),  # β1
            2π * transformtheta(initvalues3[3]) - π,  # ϕ
            initvalues3[4],  # ψt
            transformtheta(initvalues3[5]),  # βreduction1
            transformtheta(initvalues3[6]),  # βreduction2
            transformtheta(initvalues3[7])  # detection
        ],
        [
            initvalues4[1],  # R0
            transformtheta(initvalues4[2]),  # β1
            2π * transformtheta(initvalues4[3]) - π,  # ϕ
            initvalues4[4],  # ψt
            transformtheta(initvalues4[5]),  # βreduction1
            transformtheta(initvalues4[6]),  # βreduction2
            transformtheta(initvalues4[7])  # detection
        ]
    ],
)
=#
chaindict = Dict(
    "chain" => chain,
    "optimvalues" => reoptimvalues,
)

safesave(datadir("sims", "chain_omega_$(omega)_nrounds_$(n_rounds).jld2"), chaindict)




