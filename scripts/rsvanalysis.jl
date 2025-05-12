
using DrWatson

@quickactivate :ImmuneBoostingODEs

using DataFrames
using DifferentialEquations
using Memoization
using Optimization
using OptimizationOptimJL
using OptimizationPolyalgorithms
using OrdinaryDiffEq
using Random
using SciMLSensitivity
using Turing
using Zygote

#include("samplepriors.jl")

testrun = true 

n_rounds = testrun ? 25 : 10_000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function loss(
    p; 
    data, 
    prob, 
    callback, 
    saveat,
    S0=0.1,
    I0=2e-5,
    equalrs=true,
    t0=1996.737,
    alg=Vern9(; lazy=false),
    abstol=1e-15,
    maxiters=1e8,
) 
    u0 = sirns_u0_transformedp(S0, I0; p, equalrs, t0)
    sol = solve(
        prob, alg; 
        p, u0, callback, saveat, abstol, maxiters,
    )

    if !SciMLBase.successful_retcode(sol)
        @warn "$(sol.retcode) with p=$p, u0=$u0"
        return Inf
    end

    cumulativecases = modelcompartments(sol, 8)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000
    loss = sum(abs2, data .- incidentcases .* ImmuneBoostingODEs._logistic(p[7] - 4.185))
    return loss
end

function optimizesirns(
    data,
    p;
    callback,
    saveat,
    S0=0.1,
    I0=2e-5,
    equalrs=true,
    t0=1996.737,
    tspan=( 1996.737, last(saveat) ),
    alg=Vern9(; lazy=false),
    odesolverabstol=1e-15,
    odesolvermaxiters=5e7,
    optimizationsolvermaxiters=1e5,
    adtype=Optimization.AutoZygote(),
)
    u0 = sirns_u0_transformedp(S0, I0; p, equalrs, t0)
    prob = ODEProblem(transformedsirns!, u0, tspan, p)
    sol = solve(
        prob, alg; 
        p, u0, callback, saveat, abstol=odesolverabstol, maxiters=odesolvermaxiters,
    )
    
    if !SciMLBase.successful_retcode(sol)
        @error "$(sol.retcode) with p=$p, u0=$u0"
        return nothing
    end
    
    optf = Optimization.OptimizationFunction(
        (x, p) -> loss(
            x; 
            callback, 
            data, 
            prob, 
            saveat, 
            S0, 
            I0, 
            equalrs, 
            t0, 
            alg, 
            abstol=odesolverabstol, 
            maxiters=odesolvermaxiters,
        ), 
        adtype
    )
    optprob = Optimization.OptimizationProblem(optf, p)
    result_ode = Optimization.solve(optprob, PolyOpt(); maxiters=optimizationsolvermaxiters)
    return result_ode
end

initial_params1 = optimizesirns(
    data.Cases, [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ];
    callback=cbs, saveat
)

initial_params2 = optimizesirns(
    omegaspecifictransformedsirns!, data.Cases, [ 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ];
    callback=cbs, omega, saveat
)
initial_params3 = optimizesirns(
    omegaspecifictransformedsirns!, data.Cases, [ -3.0, 0.0, 0.0, 2.3, 0.0, 0.0, 0.0 ];
    callback=cbs, omega, saveat
)

tspan = ( 1996.737, last(saveat) )

initial_params1 = let
    p = [ 
        0.693,  # logr0 
        -2.20,  # logitβ1 
        0.0,  # logitϕ 
        0.0,  # logψ 
        4.60,  # logitreduction1 
        4.60,  # logitreduction2 
        -3.89,  # logitdetection
    ]
    u0 = sirns_u0_transformedp(0.01, 2e-5; p, omega, equalrs=true, t0=1996.737)
    prob = ODEProblem(omegaspecifictransformedsirns!, u0, tspan, p)
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
    )
    
    optf = Optimization.OptimizationFunction(
        (x, p) -> loss(x; callback=cbs, data=data.Cases, prob, saveat), 
        Optimization.AutoZygote()
    )
    optprob = Optimization.OptimizationProblem(optf, p)
    Optimization.solve(optprob, PolyOpt(), maxiters = 100)
end

p = [ 
    0.693,  # logr0 
    -2.20,  # logitβ1 
    0.0,  # logitϕ 
    0.0,  # logψ 
    0.0,  # logitreduction1 
    0.0,  # logitreduction2 
    0.0,  # logitdetection
]
u0 = sirns_u0_transformedp(0.01, 2e-5; p, omega, equalrs=true, t0=1996.737)
prob = ODEProblem(omegaspecifictransformedsirns!, u0, tspan, p)
sol = solve(
    prob, Vern9(; lazy=false); 
    p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
)
optf = Optimization.OptimizationFunction(
    (x, p) -> loss(x; callback=cbs, data=data.Cases, prob, saveat), 
    Optimization.AutoZygote()
)
optprob = Optimization.OptimizationProblem(optf, p)
result_ode = Optimization.solve(optprob, PolyOpt(), maxiters = 100)

initial_params1 = let
    p = [ 
        0.0,  # logr0 
        0.0,  # logitβ1 
        0.0,  # logitϕ 
        0.0,  # logψ 
        0.0,  # logitreduction1 
        0.0,  # logitreduction2 
        0.0,  # logitdetection
    ]
    u0 = sirns_u0_transformedp(0.01, 2e-5; p, omega, equalrs=true, t0=1996.737)
    prob = ODEProblem(omegaspecifictransformedsirns!, u0, tspan, p)
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
    )
    
    optf = Optimization.OptimizationFunction(
        (x, p) -> loss(x; callback=cbs, data=data.Cases, prob, saveat), 
        Optimization.AutoZygote()
    )
    optprob = Optimization.OptimizationProblem(optf, p)
    Optimization.solve(optprob, PolyOpt(), maxiters = 100)
end

initial_params2 = let
    p = [ 
        3.0,  # logr0 
        0.0,  # logitβ1 
        0.0,  # logitϕ 
        0.0,  # logψ 
        0.0,  # logitreduction1 
        0.0,  # logitreduction2 
        0.0,  # logitdetection
    ]
    u0 = sirns_u0_transformedp(0.01, 2e-5; p, omega, equalrs=true, t0=1996.737)
    prob = ODEProblem(omegaspecifictransformedsirns!, u0, tspan, p)
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
    )
    optf = Optimization.OptimizationFunction(
        (x, p) -> loss(x; callback=cbs, data=data.Cases, prob, saveat), 
        Optimization.AutoZygote()
    )
    optprob = Optimization.OptimizationProblem(optf, p)
    result_ode = Optimization.solve(optprob, PolyOpt(); maxiters=100)
    result_ode
end


# example parameters and initial conditions 
p = [ 
    0.693,  # logr0 
    -2.20,  # logitβ1 
    0.0,  # logitϕ 
    0.0,  # logψ 
    4.60,  # logitreduction1 
    4.60,  # logitreduction2 
    -3.89,  # logitdetection
]
u0 = sirns_u0_transformedp(0.01, 2e-5; p, omega=2.0, equalrs=true, t0=1996.737)


prob = ODEProblem(transformedsirns_omega2!, u0, tspan, p)
sol = solve(
    prob, Vern9(; lazy=false); 
    p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
)



adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction(
    (x, p) -> loss(x; callback=cbs, data=data.Cases, prob, saveat), 
    adtype
)
optprob = Optimization.OptimizationProblem(optf, p)

result_ode = Optimization.solve(optprob, PolyOpt(), maxiters = 100)

initial_params = let
    gamma=48.7; mu=0.0087
    logr0, logitβ1, logitϕ, logψ, logitreduction1, logitreduction2, logitdetection, = result_ode
    rzero = exp(logr0)
    betazero = rzero * (gamma + mu)
    betaone = ImmuneBoostingODEs._logistic(logitβ1)
    phi = ImmuneBoostingODEs._logistic(logitϕ) * 2π - π
    psi_t = logψ - log(0.7)
    reduction1 = ImmuneBoostingODEs._logistic(logitreduction1)
    reduction2 = ImmuneBoostingODEs._logistic(logitreduction2)
    detection = ImmuneBoostingODEs._logistic(logitdetection)
    [
        [ 
            rzero, 
            betazero, 
            phi, 
            psi_t, 
            reduction1, 
            reduction2, 
            [ 0.05, 0.1, 0.2, 0.5 ][i],  # pparameter_t
            detection 
        ]
        for i ∈ 1:4
    ]
    
end

prob = fittedsimulationsetup(saveat)

include("rsvfitmodel.jl")

chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat; omega),
    Turing.NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    #16;
    #initial_params=optimvalues,
    initial_params,
)

chaindf = DataFrame(chain)
plotchains(chaindf)

chaindict = Dict(
    "chain" => chain,
    "optimvalues" => reoptimvalues,
)

safesave(datadir("sims", "chain_omega_$(omega)_nrounds_$(n_rounds).jld2"), chaindict)
