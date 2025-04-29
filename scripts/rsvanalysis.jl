
using DrWatson

@quickactivate :ImmuneBoostingODEs
using DataFrames, DifferentialEquations, Optim, Random, Turing

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

transformtheta(x) = exp(x) / (1 + exp(x))

function optimmodel(paramvector, β1, ω, data; callback, saveat)
    # take any real values in paramvector and transform into "acceptable" values 
    R0, ψ, _detection = paramvector
    detection = transformtheta(_detection)
    if R0 < 0 || ψ < 0 || detection < 0 || detection > 1
        return Inf
    end

    β0 = R0 * 48.7087
    p = SirnsParameters(
        β0, 
        β1, 
        0.0,  # ϕ 
        48.7,  # γ 
        0.0087,  # μ 
        ψ, 
        ω, 
        β0, 
        0.8 * β0, 
        0.9 * β0
    )
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection
    sol = memosolver(
        prob, Vern9(; lazy=false); 
        p, u0, callback, saveat, save_idxs=[ 8 ], 
        abstol=1e-15, maxiters=1e8, verbose=false,
    )
    if sol.retcode != :Success
        return Inf
    end
    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000 .* detection
    return sum([ abs2(data[i] - incidentcases[i]) for i ∈ eachindex(data) ])
end

result0001 = optimize(
    x -> optimmodel(x, 0.001, omega, data.Cases; callback=cbs, saveat), 
    [ 2.0, 0.1, -4 ], 
    NelderMead(), 
    Optim.Options(; iterations=50_000)
)
result001 = optimize(
    x -> optimmodel(x, 0.01, omega, data.Cases; callback=cbs, saveat), 
    [ 2.0, 0.1, -4 ], 
    NelderMead(), 
    Optim.Options(; iterations=50_000)
)
result01 = optimize(
    x -> optimmodel(x, 0.1, omega, data.Cases; callback=cbs, saveat), 
    [ 2.0, 0.1, -4 ], 
    NelderMead(), 
    Optim.Options(; iterations=50_000)
)
result025 = optimize(
    x -> optimmodel(x, 0.25, omega, data.Cases; callback=cbs, saveat), 
    [ 2.0, 0.1, -4 ], 
    NelderMead(), 
    Optim.Options(; iterations=50_000)
)

initvalues0001 = Optim.minimizer(result0001)
initvalues001 = Optim.minimizer(result001)
initvalues01 = Optim.minimizer(result01)
initvalues025 = Optim.minimizer(result025)

Random.seed!(round(Int, omega * 100))

chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat; omega),
    NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    initial_params=[
        [
            initvalues0001[1],  # R0
            0.001,  # β1
            0.0,  # ϕ
            initvalues0001[2],  # ψ
            0.8,  # βreduction1
            0.9,  # βreduction2
            transformtheta(initvalues0001[3])
        ],
        [
            initvalues001[1],  # R0
            0.01,  # β1
            0.0,  # ϕ
            initvalues001[2],  # ψ
            0.8,  # βreduction1
            0.9,  # βreduction2
            transformtheta(initvalues001[3])
        ],
        [
            initvalues01[1],  # R0
            0.1,  # β1
            0.0,  # ϕ
            initvalues01[2],  # ψ
            0.8,  # βreduction1
            0.9,  # βreduction2
            transformtheta(initvalues01[3])
        ],
        [
            initvalues025[1],  # R0
            0.25,  # β1
            0.0,  # ϕ
            initvalues025[2],  # ψ
            0.8,  # βreduction1
            0.9,  # βreduction2
            transformtheta(initvalues025[3])
        ]
    ],
)

chaindict = Dict(
    "chain" => chain,
    "result0001" => result0001,
    "result001" => result001,
    "result01" => result01,
    "result025" => result025,
)

safesave(datadir("sims", "chain_omega_$(omega)_nrounds_$(n_rounds).jld2"), chaindict)
