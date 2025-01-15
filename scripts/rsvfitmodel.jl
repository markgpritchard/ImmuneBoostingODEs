
@model function fitmodel(
    incidence, prob, cbs, saveat;
    betazeroprior=truncated(Exponential(150), 0, 4870),  # truncated at R0 = 100
    betaoneprior=Uniform(0, 0.9),
    phiprior=Uniform(-π, π),
    psiprior=truncated(Exponential(1), 0, 1000),
    betareduction1prior=Beta(4, 1),
    betareduction2prior=Beta(9, 1),
    omega=2.0,
    detectionprior=Beta(1, 98)
)
    β0 ~ betazeroprior
    β1 ~ betaoneprior
    ϕ ~ phiprior
    γ = 48.7  # generation time 7.5 days
    μ = 0.0087  # Scotland's birth rate = 48000 / 5.5e6
    ψ ~ psiprior
    ω = omega
    βreduction1 ~ betareduction1prior
    βreduction2 ~ betareduction2prior
    detection ~ detectionprior

    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, β0, βreduction1 * β0, βreduction2 * β0)
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection

    sol = memosolver(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, save_idxs=[ 8 ], 
        abstol=1e-15, maxiters=1e8, verbose=false,
    )
    if sol.retcode != :Success
        #@info "Adding logprob -Inf when p=$p, detection=$detection"
        Turing.@addlogprob! -Inf
        return nothing
    end

    #cumulativecases = modelcompartments(sol, :cc)
    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000 .* detection

    for i ∈ eachindex(incidentcases)
        incidence[i] ~ Poisson(incidentcases[i] + 1e-10)
    end
end
