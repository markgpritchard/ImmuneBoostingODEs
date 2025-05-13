
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

testrun = false 

n_rounds = testrun ? 25 : 10_000
optimizationsolvermaxiters = testrun ? 100 : 1e5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting parameters 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# `parms` is a vector containing the following in order:
    # r0, 
    # logitβ1, 
    # ϕ, 
    # logψ, 
    # logω, 
    # logitbetaprimemultiplier 
    # logitfinalbetaprime, 
    # logitproportiondetected, 
    # and finally any number, which will be mutated during Optimization (not related to
        # `rparameter` in `fitmodel`) 
function loss(
    parms; 
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
    gamma=48.7,
    mu=0.0087,
) 
    if parms[1] < 0 || parms[3] < -π || parms[3] > π
        return Inf
    end

    u0 = sirns_u0_transformedp(S0, I0; p=parms, equalrs, t0)
    sol = solve(
        prob, alg; 
        p=parms, u0, callback, saveat, abstol, maxiters,
    )

    if !SciMLBase.successful_retcode(sol)
        @warn "$(sol.retcode) with parms=$parms, u0=$u0"
        return Inf
    end

    cumulativecases = modelcompartments(sol, 8)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000
    loss = sum(abs2, data .- incidentcases .* ImmuneBoostingODEs._logistic(parms[7]))
    #loss = sum(abs2, data .- incidentcases .* 0.015)  # fix at 1.5% detected
    #println("parms=$parms -> totalcases=$(last(cumulativecases)), loss=$loss")
    return loss
end

function optimizesirns(
    data,
    parms;
    callback,
    saveat,
    S0=0.1,
    I0=2e-5,
    equalrs=true,
    t0=1996.737,
    tspan=( 1996.737, last(saveat) ),
    gamma=48.7,
    mu=0.0087,
    alg=Vern9(; lazy=false),
    odesolverabstol=1e-15,
    odesolvermaxiters=5e7,
    optimizationsolvermaxiters=1e5,
    adtype=Optimization.AutoZygote(),
)
    u0 = sirns_u0_transformedp(S0, I0; p=parms, equalrs, t0)
    prob = ODEProblem(transformedsirns!, u0, tspan, parms)
    sol = solve(
        prob, alg; 
        p=parms, u0, callback, saveat, abstol=odesolverabstol, maxiters=odesolvermaxiters,
    )
    
    if !SciMLBase.successful_retcode(sol)
        @error "$(sol.retcode) with parms=$parms, u0=$u0"
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
    optprob = Optimization.OptimizationProblem(optf, parms)
    result_ode = Optimization.solve(optprob, PolyOpt(); maxiters=optimizationsolvermaxiters)
    return result_ode
end

initial_params1 = optimizesirns(
    data.Cases, 
    [ 2.0, log(0.1), 0.0, log(1), log(2), log(1), log(0.9 / 0.1), log(0.02), 1.0 ];
    callback=optimcbs, saveat, optimizationsolvermaxiters,
)
initial_params2 = optimizesirns(
    data.Cases, 
    [ 1.0, log(1), -0.5, log(0.001), log(1), log(0.1), log(0.8 / 0.2), log(0.015), 1.0 ];
    callback=optimcbs, saveat, optimizationsolvermaxiters,
)
initial_params3 = optimizesirns(
    data.Cases, 
    [ 5.0, log(0.2), 0.5, log(10), log(0.5), log(0.99/0.01), log(0.9/0.1), log(0.01), 1.0 ];
    callback=optimcbs, saveat, optimizationsolvermaxiters,
)
initial_params4 = optimizesirns(
    data.Cases, 
    [ 10.0, log(0.001), 0.0, log(0.5), log(2), log(1), log(0.99 / 0.01), log(0.005), 1.0 ];
    callback=optimcbs, saveat, optimizationsolvermaxiters,
)

tspan = ( 1996.737, last(saveat) )
initialp = SirnsParameters(
    2.0 * (48.7 + 0.0087),  # β0::T
    0.1,  # β1::T
    0.0,  # ϕ::T
    48.7,  # γ::Float64
    0.0087,  # μ::Float64 
    1.0,  # ψ::T
    1.0,  # ω::T
    2.0 * (48.7 + 0.0087),  # originalβ0::T
    0.5,  # betaprimemultiplier::T
    0.5,  # finalbetaprime::T
    0.5,  # proportiondetected::T
)
u0 = sirns_u0(0.01, 2e-5; p=initialp, equalrs=true, t0=1996.737)
prob = ODEProblem(sirns!, u0, tspan, initialp)

include("rsvfitmodel.jl")

chain = sample(
    fitmodel(data.Cases, prob, cbs, saveat),
    Turing.NUTS(0.65),
    MCMCThreads(),
    n_rounds,
    4;
    #16;
    #initial_params=optimvalues,
    initial_params=[
        [ initial_params1[1:8]; 1.0 ],
        [ initial_params2[1:8]; 1.0 ],
        [ initial_params3[1:8]; 1.0 ],
        [ initial_params4[1:8]; 1.0 ],
    ]
)

chaindf = DataFrame(chain)
plotchains(chaindf)

chaindict = Dict(
    "chain" => chain,
    "initial_params1" => initial_params1,
    "initial_params2" => initial_params2,
    "initial_params3" => initial_params3,
    "initial_params4" => initial_params4,
    "n_rounds" => n_rounds,
    "optimizationsolvermaxiters" => optimizationsolvermaxiters,
)

safesave(datadir("sims", "chaindict_nrounds_$(n_rounds).jld2"), chaindict)


##########
#=
modeloutputs = Vector{Vector{Float64}}(undef, size(chaindf, 1))

for i ∈ axes(chaindf, 1)
    p = SirnsParameters(
        chaindf.r0[i] / (48.7 + 0.0087),  # β0::T
        _logistic(chaindf.logitβ1[i]),  # β1::T
        chaindf.ϕ[i],  # ϕ::T
        48.7,  # γ::Float64
        0.0087,  # μ::Float64 
        exp(chaindf.logψ[i]),  # ψ::T
        exp(chaindf.logω[i]),  # ω::T
        chaindf.r0[i] / (48.7 + 0.0087),  # originalβ0::T
        _logistic(chaindf.logitfinalbetaprime[i]),  # finalbetaprime::T
        _logistic(chaindf.logitproportiondetected[i]),  # proportiondetected::T
    )
    u0 = sirns_u0(0.1, 2e-5; p, t0=1996.737)
    tspan = ( 1996.737, last(saveat) )
    
    prob = ODEProblem(sirns!, u0, tspan, p) 
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8,
    )
    cumulativecases = modelcompartments(sol, 8)
    incidentcases = casespertimeblock(cumulativecases .* 5_450_000 * p.proportiondetected)
    modeloutputs[i] = incidentcases
end

medianoutput = Vector{Float64}(undef, 353)
lcioutput = Vector{Float64}(undef, 353)
ucioutput = Vector{Float64}(undef, 353)
for t ∈ 1:353 
    lc, me, uc = quantile([ modeloutputs[i][t] for i ∈ axes(chaindf, 1) ], [ 0.025, 0.5, 0.975 ])
    medianoutput[t] = me
    lcioutput[t] = lc 
    ucioutput[t] = uc 
end


fig = Figure(; size=( 500, 500 ))

ga = GridLayout(fig[1, 1])
ax = Axis(ga[1, 1]; xticks=2017:2:2023, yticks=0:200:600)
lines!(ax, data.Date, medianoutput; color=COLOURVECTOR[1], linewidth=1,)
band!(ax, data.Date, lcioutput, ucioutput; color=( COLOURVECTOR[1], 0.5 ))
scatter!(ax, data.Date, data.Cases; color=:black, markersize=3)


fig








for (i, v) ∈ enumerate(plotvvector)
    plotfittedsimulationquantiles!(axs[i], data, v, saveat)
    text!(
        axs[i], textlocation[1], textlocation[2]; 
        text="ω=$(omegalabels[i])", fontsize=11.84, align=( :left, :top )
    )
end

stringencyax = Axis(ga[1:7, 1])
vspan!(stringencyax, reduceday, increaseday, color=( :gray, 0.1 ))
for x ∈ 2017:1:2023
    vlines!(
        stringencyax, x; 
        color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
    )
end

gb = GridLayout(fig[1, 2])
ax2 = Axis(
    gb[1, 1]; 
    xticks=( logomegavalues, omegalabels ), 
    yticks=( log.([ 1, 2, 5, 10, 20, 40 ]), [ "1", "2", "5", "10", "20", "40" ])
)
scatter!(
    ax2, 
    logomegavalues, 
    log.([ quantile(v.β0, 0.5) for v ∈ pv ] ./ (γ + μ)); 
    color=:blue, markersize=5,
)
rangebars!(
    ax2, 
    logomegavalues, 
    log.([ quantile(v.β0, 0.05) for v ∈ pv ] ./ (γ + μ)), 
    log.([ quantile(v.β0, 0.95) for v ∈ pv ] ./ (γ + μ));
    color=:blue,
)
for y ∈ [ 1, 2, 5, 10, 20, 40 ]
    hlines!(
        ax2, log(y); 
        color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
    )
end
ax3 = Axis(gb[2, 1]; xticks=( logomegavalues, omegalabels ), yticks=0:5:20,)
scatter!(
    ax3, 
    logomegavalues, 
    #[ quantile(v.β1, 0.5) for v ∈ pv ] .* [ quantile(v.β0, 0.5) for v ∈ pv ] ./ (γ + μ); 
    100 .* [ quantile(v.β1, 0.5) for v ∈ pv ]; 
    color=:blue, markersize=5,
)
for y ∈ 0:5:20
    hlines!(
        ax3, y; 
        color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
    )
end    
rangebars!(
    ax3, 
    logomegavalues, 
    #[ quantile(v.β1, 0.05) for v ∈ pv ] .* [ quantile(v.β0, 0.05) for v ∈ pv ] ./ (γ + μ), 
    #[ quantile(v.β1, 0.95) for v ∈ pv ] .* [ quantile(v.β0, 0.95) for v ∈ pv ] ./ (γ + μ);
    100 .* [ quantile(v.β1, 0.05) for v ∈ pv ], 
    100 .* [ quantile(v.β1, 0.95) for v ∈ pv ];
    color=:blue,
)
ax4 = Axis(
    gb[3, 1]; 
    xticks=( logomegavalues, omegalabels ), 
    yticks=( 
        log.([ 0.001, 0.1, 10, 1000 ]), 
        [ "0.001", "0.1", "10", "1000" ]
    )
)
scatter!(
    ax4, logomegavalues, log.([ quantile(v.ψ, 0.5) for v ∈ pv ]); 
    color=:blue, markersize=5,
)
rangebars!(
    ax4, 
    logomegavalues, 
    log.([ quantile(v.ψ, 0.05) for v ∈ pv ]), 
    log.([ quantile(v.ψ, 0.95) for v ∈ pv ]);
    color=:blue,
)
for y ∈ [ 0.001, 0.1, 10, 1000 ]
    hlines!(
        ax4, log(y); 
        color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
    )
end   
ax5 = Axis(gb[4, 1]; xticks=( logomegavalues, omegalabels ), yticks=20:10:50)
scatter!(
    ax5, logomegavalues, 100 .* (1 .- [ quantile(v.βreduction1, 0.5) for v ∈ pv ]); 
    color=:blue, markersize=5,
)
rangebars!(
    ax5, 
    logomegavalues, 
    100 .* (1 .- [ quantile(v.βreduction1, 0.05) for v ∈ pv ]), 
    100 .* (1 .- [ quantile(v.βreduction1, 0.95) for v ∈ pv ]);
    color=:blue,
)
for y ∈ 20:10:50
    hlines!(
        ax5, y; 
        color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
    )
end   
ax6 = Axis(gb[5, 1]; xticks=( logomegavalues, omegalabels ), yticks=0.0:0.5:2.0)
scatter!(
    ax6, logomegavalues, [ quantile(v.detection, 0.5) for v ∈ pv ] .* 100; 
    color=:blue, markersize=5,
)
rangebars!(
    ax6, 
    logomegavalues, 
    [ quantile(v.detection, 0.05) for v ∈ pv ] .* 100, 
    [ quantile(v.detection, 0.95) for v ∈ pv ] .* 100;
    color=:blue,
)
for y ∈ 0.0:0.5:2.0
    hlines!(
        ax6, y; 
        color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
    )
end   

linkxaxes!(stringencyax, axs...)
for i ∈ 1:7 
    formataxis!(
        axs[i]; 
        hidex=(i != 7), hidexticks=(i != 7), trimspines=true, hidespines=( :t, :r ),
        setpoint=textlocation,
    )
    if i != 7 hidespines!(axs[i], :b) end
end
formataxis!(
    stringencyax; 
    hidespines=( :l, :r, :t, :b ), 
    hidex=true, hidexticks=true, hidey=true, hideyticks=true
)
Label(
    ga[1:7, 0], "Weekly incidence"; 
    fontsize=11.84, rotation=π/2, tellheight=false
)
Label(ga[8, 1], "Year"; fontsize=11.84, tellwidth=false)
colgap!(ga, 1, 5)
for r ∈ [ 1, 9 ] rowgap!(ga, 7, 5) end
for (i, ax) ∈ enumerate([ ax2, ax3, ax4, ax5, ax6 ])
    formataxis!(
        ax; 
        hidex=(i != 5), hidexticks=(i != 5), trimspines=true, hidespines=( :t, :r ),
    )
    if i != 5 hidespines!(ax, :b) end
    if i != 4 setvalue!(ax, 1, 0) end
end
Label(gb[1, 0], L"Mean $\mathcal{R}_0$"; fontsize=11.84, rotation=π/2, tellheight=false)
Label(
    gb[2, 0], "Magnitude of\n forcing, %"; 
    fontsize=11.84, rotation=π/2, tellheight=false
)
Label(
    gb[3, 0], L"$\psi$"; 
    fontsize=11.84, rotation=π/2, tellheight=false
)
Label(
    gb[4, 0], "Effect of\ninterventions, %"; 
    fontsize=11.84, rotation=π/2, tellheight=false
)
Label(
    gb[5, 0], "Proportion\ndiagnosed, %"; 
    fontsize=11.84, rotation=π/2, tellheight=false
)
Label(gb[6, 1], "Waning rate, ω"; fontsize=11.84, tellwidth=false)
colgap!(gb, 1, 5)
rowgap!(gb, 5, 5)

labelplots!([ "A", "B", ], [ ga, gb ]; rows=[ 1, 1 ])

fig
=#
