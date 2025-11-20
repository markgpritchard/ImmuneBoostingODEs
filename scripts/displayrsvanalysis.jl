
using DrWatson
#@quickactivate "ImmuneBoostingODEs"
@quickactivate :ImmuneBoostingODEs

using CairoMakie
using DataFrames
using DifferentialEquations
using Random
using Turing


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("rsvsetup.jl")

println("In the 12 months from 1 April each year")
for y ∈ 2016:2022 
    inds = findall(x -> y <= x < y + 1, data.AprilYear)
    println("    $(sum(data.Cases[inds])) cases in $y")
end 

prob = fittedsimulationsetup(saveat)

model_psi0 = fitmodel(data.Cases, prob; callback=betareductioncallback, saveat, psi=0)
priors_psi0 = sample(model_psi0, Prior(), MCMCThreads(), 1_000, 4)

predmodel = fitmodel(missing, prob; callback=betareductioncallback, saveat, psi=0)
predictions = predict(Random.default_rng(), predmodel, priors_psi0)
predictionarray = Array(predictions)
predictionquantiles = zeros(length(data.Cases), 7)
size(predictionarray, 2) == size(predictionquantiles, 1)
for i in axes(predictionquantiles, 1)
    predictionquantiles[i, :] .= max.(
        0, 
        quantile(
            skipmissing(predictionarray[:, i]), [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
        )
    )
end

fig = let 
    fig = Figure()
    axs = Axis(fig[1, 1])
    band!(data.Date, predictionquantiles[:, 1], predictionquantiles[:, 7]; color=(COLOUR_I, 0.3),)
    band!(data.Date, predictionquantiles[:, 2], predictionquantiles[:, 6]; color=(COLOUR_I, 0.5),)
    band!(data.Date, predictionquantiles[:, 3], predictionquantiles[:, 5]; color=(COLOUR_I, 0.7),)
    lines!(data.Date, predictionquantiles[:, 4]; color=COLOUR_I, linewidth=1,)
    scatter!(data.Date, data.Cases; color=:black, marker=:x, markersize=3,)

    fig
end

initial_params_psi0 = initialparams_map(model_psi0, 4; maxiters=10,)
#=
paramnames = (:β0, :β1, :ϕ, :ω, :βreductionfactor, :detection, :minsigma2, :S0, :I0, :R1, :R2)
Threads.@threads for k in 1:4 
    ip = [rand(Uniform(lowerbounds[i], upperbounds[i])) for i in eachindex(lowerbounds)]
    ipf = optimizesirns(
        data.Cases, prob, ip;
        callback=betareductioncallback,
        saveat,  
        lb=lowerbounds, 
        ub=upperbounds, 
        verbosity=0, 
        psi=0,
        optimizationsolvermaxiters=5000,  # to increase later
    )
    #initial_params_psi0[k] = exp.(ipf.minimizer)
    paramvalues = Tuple(exp.(ipf.minimizer))
    if paramvalues[3] > π 
        paramvalues[3] += -2π 
    end
    namedtup = NamedTuple{paramnames}(paramvalues)
    initial_params_psi0[k] = InitFromParams(namedtup)
    @info "initial_params[$k] ($(ipf.retcode)) = $(initial_params_psi0[k])"
end

=#

#=
# which is a vector containing the log of the following parameters in order:
        # β0 
        # β1 
        # ϕ 
        # ω 
        # βreductionfactor 
        # detection 
        # minsigma2 
        # S0 
        # I0 
        # R1 
        # R2
        =#

#maximum_a_posteriori(model_psi0)

samples_psi0 = sample(model_psi0, NUTS(100, 0.65), MCMCThreads(), 100, 4; initial_params=initial_params_psi0)

predmodel = fitmodel(missing, prob; callback=betareductioncallback, saveat, psi=0)
predictions = predict(Random.default_rng(), predmodel, samples_psi0)
predictionarray = Array(predictions)
predictionquantiles = zeros(length(data.Cases), 7)
size(predictionarray, 2) == size(predictionquantiles, 1)
for i in axes(predictionquantiles, 1)
    predictionquantiles[i, :] .= max.(0, quantile(skipmissing(predictionarray[:, i]), [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]))
end

fig = let 
    fig = Figure()
    axs = Axis(fig[1, 1])
    band!(data.Date, predictionquantiles[:, 1], predictionquantiles[:, 7]; color=(COLOUR_I, 0.3),)
    band!(data.Date, predictionquantiles[:, 2], predictionquantiles[:, 6]; color=(COLOUR_I, 0.5),)
    band!(data.Date, predictionquantiles[:, 3], predictionquantiles[:, 5]; color=(COLOUR_I, 0.7),)
    lines!(data.Date, predictionquantiles[:, 4]; color=COLOUR_I, linewidth=1,)
    scatter!(data.Date, data.Cases; color=:black, marker=:x, markersize=3,)

    fig
end



#

model_psi05 = fitmodel(data.Cases, prob; callback=betareductioncallback, saveat, psi=0.5)
priors_psi05 = sample(model_psi05, Prior(), MCMCThreads(), 1_000, 4)

predmodel = fitmodel(missing, prob; callback=betareductioncallback, saveat, psi=0.5)
predictions = predict(Random.default_rng(), predmodel, priors_psi05)
predictionarray = Array(predictions)
predictionquantiles = zeros(length(data.Cases), 7)
size(predictionarray, 2) == size(predictionquantiles, 1)
for i in axes(predictionquantiles, 1)
    predictionquantiles[i, :] .= quantile(skipmissing(predictionarray[:, i]), [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975])
end

fig = let 
    fig = Figure()
    axs = Axis(fig[1, 1])
    band!(data.Date, predictionquantiles[:, 1], predictionquantiles[:, 7]; color=(COLOUR_I, 0.3),)
    band!(data.Date, predictionquantiles[:, 2], predictionquantiles[:, 6]; color=(COLOUR_I, 0.5),)
    band!(data.Date, predictionquantiles[:, 3], predictionquantiles[:, 5]; color=(COLOUR_I, 0.7),)
    lines!(data.Date, predictionquantiles[:, 4]; color=COLOUR_I, linewidth=1,)
    scatter!(data.Date, data.Cases; color=:black, marker=:x, markersize=3,)

    fig
end

model_psi1 = fitmodel(data.Cases, prob; callback=betareductioncallback, saveat, psi=1)
priors_psi1 = sample(model_psi1, Prior(), MCMCThreads(), 1_000, 4)

predmodel = fitmodel(missing, prob; callback=betareductioncallback, saveat, psi=1)
predictions = predict(Random.default_rng(), predmodel, priors_psi1)
predictionarray = Array(predictions)
predictionquantiles = zeros(length(data.Cases), 7)
size(predictionarray, 2) == size(predictionquantiles, 1)
for i in axes(predictionquantiles, 1)
    predictionquantiles[i, :] .= quantile(skipmissing(predictionarray[:, i]), [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975])
end

fig = let 
    fig = Figure()
    axs = Axis(fig[1, 1])
    band!(data.Date, predictionquantiles[:, 1], predictionquantiles[:, 7]; color=(COLOUR_I, 0.3),)
    band!(data.Date, predictionquantiles[:, 2], predictionquantiles[:, 6]; color=(COLOUR_I, 0.5),)
    band!(data.Date, predictionquantiles[:, 3], predictionquantiles[:, 5]; color=(COLOUR_I, 0.7),)
    lines!(data.Date, predictionquantiles[:, 4]; color=COLOUR_I, linewidth=1,)
    scatter!(data.Date, data.Cases; color=:black, marker=:x, markersize=3,)

    fig
end

model_psi5 = fitmodel(data.Cases, prob; callback=betareductioncallback, saveat, psi=5)
priors_psi5 = sample(model_psi5, Prior(), MCMCThreads(), 1_000, 4)

predmodel = fitmodel(missing, prob; callback=betareductioncallback, saveat, psi=5)
predictions = predict(Random.default_rng(), predmodel, priors_psi5)
predictionarray = Array(predictions)
predictionquantiles = zeros(length(data.Cases), 7)
size(predictionarray, 2) == size(predictionquantiles, 1)
for i in axes(predictionquantiles, 1)
    predictionquantiles[i, :] .= quantile(skipmissing(predictionarray[:, i]), [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975])
end

fig = let 
    fig = Figure()
    axs = Axis(fig[1, 1])
    band!(data.Date, predictionquantiles[:, 1], predictionquantiles[:, 7]; color=(COLOUR_I, 0.3),)
    band!(data.Date, predictionquantiles[:, 2], predictionquantiles[:, 6]; color=(COLOUR_I, 0.5),)
    band!(data.Date, predictionquantiles[:, 3], predictionquantiles[:, 5]; color=(COLOUR_I, 0.7),)
    lines!(data.Date, predictionquantiles[:, 4]; color=COLOUR_I, linewidth=1,)
    scatter!(data.Date, data.Cases; color=:black, marker=:x, markersize=3,)

    fig
end

model_psi10 = fitmodel(data.Cases, prob; callback=betareductioncallback, saveat, psi=10)
priors_psi10 = sample(model_psi10, Prior(), MCMCThreads(), 1_000, 4)

predmodel = fitmodel(missing, prob; callback=betareductioncallback, saveat, psi=10)
predictions = predict(Random.default_rng(), predmodel, priors_psi10)
predictionarray = Array(predictions)
predictionquantiles = zeros(length(data.Cases), 7)
size(predictionarray, 2) == size(predictionquantiles, 1)
for i in axes(predictionquantiles, 1)
    predictionquantiles[i, :] .= quantile(skipmissing(predictionarray[:, i]), [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975])
end

fig = let 
    fig = Figure()
    axs = Axis(fig[1, 1])
    band!(data.Date, predictionquantiles[:, 1], predictionquantiles[:, 7]; color=(COLOUR_I, 0.3),)
    band!(data.Date, predictionquantiles[:, 2], predictionquantiles[:, 6]; color=(COLOUR_I, 0.5),)
    band!(data.Date, predictionquantiles[:, 3], predictionquantiles[:, 5]; color=(COLOUR_I, 0.7),)
    lines!(data.Date, predictionquantiles[:, 4]; color=COLOUR_I, linewidth=1,)
    scatter!(data.Date, data.Cases; color=:black, marker=:x, markersize=3,)

    fig
end

model_psi20 = fitmodel(data.Cases, prob; callback=betareductioncallback, saveat, psi=20)
priors_psi20 = sample(model_psi20, Prior(), MCMCThreads(), 1_000, 4)

predmodel = fitmodel(missing, prob; callback=betareductioncallback, saveat, psi=20)
predictions = predict(Random.default_rng(), predmodel, priors_psi20)
predictionarray = Array(predictions)
predictionquantiles = zeros(length(data.Cases), 7)
size(predictionarray, 2) == size(predictionquantiles, 1)
for i in axes(predictionquantiles, 1)
    predictionquantiles[i, :] .= quantile(skipmissing(predictionarray[:, i]), [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975])
end

fig = let 
    fig = Figure()
    axs = Axis(fig[1, 1])
    band!(data.Date, predictionquantiles[:, 1], predictionquantiles[:, 7]; color=(COLOUR_I, 0.3),)
    band!(data.Date, predictionquantiles[:, 2], predictionquantiles[:, 6]; color=(COLOUR_I, 0.5),)
    band!(data.Date, predictionquantiles[:, 3], predictionquantiles[:, 5]; color=(COLOUR_I, 0.7),)
    lines!(data.Date, predictionquantiles[:, 4]; color=COLOUR_I, linewidth=1,)
    scatter!(data.Date, data.Cases; color=:black, marker=:x, markersize=3,)

    fig
end



##








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load results 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rsvparameters01 = loadrsvdata(0.1) 
plotchains(rsvparameters01)
#filter!(:chain => x -> x ∈ [ 2, 4, 5 ], rsvparameters01)
plotchains(rsvparameters01)
plotvals01 = fittedsimulationquantiles(rsvparameters01, 0.1, saveat, cbs)
R0_01 = let 
    beta0s = quantile(rsvparameters01.β0, [ 0.05, 0.5, 0.95 ])
    beta0s ./ (48.7 + 0.0087)
end
#3-element Vector{Float64}:
# 27.9398281195748
# 29.58401062939877
# 30.730225320740598

rsvparameters02 = loadrsvdata(0.2)
plotchains(rsvparameters02)
filter!(:chain => x -> x == 1, rsvparameters02)
plotchains(rsvparameters02)
plotvals02 = fittedsimulationquantiles(rsvparameters02, 0.2, saveat, cbs)

rsvparameters04 = loadrsvdata(0.4)
plotchains(rsvparameters04)
filter!(:chain => x -> x == 1, rsvparameters04)
plotchains(rsvparameters04)
plotvals04 = fittedsimulationquantiles(rsvparameters04, 0.4, saveat, cbs)
βreduction1_04 = 1 .- quantile(rsvparameters04.βreduction1, [ 0.95, 0.5, 0.05 ])
#3-element Vector{Float64}:
# 0.42537389746949894
# 0.4274937278262533
# 0.4304722062246449

rsvparameters1 = loadrsvdata(1.0)
plotchains(rsvparameters1)
filter!(:chain => x -> x != 1, rsvparameters1)
plotchains(rsvparameters1)
plotvals1 = fittedsimulationquantiles(rsvparameters1, 1.0, saveat, cbs)
beta1_1 = quantile(rsvparameters1.β1, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.1754460107922378
# 0.17788693495654073
# 0.1818164970772631
psi_1 = quantile(rsvparameters1.ψ, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.00010909782961219074
# 0.0013969680682125209
# 0.005838432275809601

rsvparameters2 = loadrsvdata(2.0)
plotchains(rsvparameters2)
filter!(:chain => x -> x ∈ [ 2, 3, 5 ], rsvparameters2)
plotchains(rsvparameters2)
plotvals2 = fittedsimulationquantiles(rsvparameters2, 2.0, saveat, cbs)
βreduction1_2 = 1 .- quantile(rsvparameters2.βreduction1, [ 0.95, 0.5, 0.05 ])
#3-element Vector{Float64}:
# 0.27889082920822683
# 0.2816906091386049
# 0.28507902779791494

rsvparameters4 = loadrsvdata(4.0)
plotchains(rsvparameters4)
filter!(:chain => x -> x ∈ [ 2, 3, 5 ], rsvparameters4)
plotchains(rsvparameters4)
plotvals4 = fittedsimulationquantiles(rsvparameters4, 4.0, saveat, cbs)

rsvparameters6 = loadrsvdata(6.0)
plotchains(rsvparameters6)
filter!(:chain => x -> x != 1, rsvparameters6)
plotchains(rsvparameters6)
plotvals6 = fittedsimulationquantiles(rsvparameters6, 6.0, saveat, cbs)
R0_6 = let 
    beta0s = quantile(rsvparameters6.β0, [ 0.05, 0.5, 0.95 ])
    beta0s ./ (48.7 + 0.0087)
end
#3-element Vector{Float64}:
# 1.5808061809771425
# 1.592694730925088
# 1.6042283249278522
beta1_6 = quantile(rsvparameters6.β1, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 0.06894139362433213
# 0.07223981399834295
# 0.07529313242053985
psi_6 = quantile(rsvparameters6.ψ, [ 0.05, 0.5, 0.95 ])
#3-element Vector{Float64}:
# 288.26257743152075
# 301.2348000922061
# 314.895983917204
