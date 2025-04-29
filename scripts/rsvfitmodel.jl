
@model function fitmodel(
    incidence, prob, cbs, saveat;
    rzeroprior=Gamma(2/3, 3),
    betaoneprior=Beta(1, 4),
    phiprior=truncated(Normal(0, π/2), -π, π),
    psiprior_t=TDist(2),
    betareduction1prior=Beta(4, 1),
    betareduction2prior=Beta(9, 1),
    omega,
    detectionprior=Beta(1, 49)
)
    mean_R0 ~ rzeroprior
    β1 ~ betaoneprior
    ϕ ~ phiprior
    ψt ~ psiprior_t
    ω = omega
    βreduction1 ~ betareduction1prior
    βreduction2 ~ betareduction2prior
    detection ~ detectionprior

    γ = 48.7  # generation time 7.5 days
    μ = 0.0087  # Scotland's annual birth rate = 48000 / 5.5e6
    β0 = mean_R0 * 48.7087
    ψ = exp(0.7 * ψt) 
    p = SirnsParameters(β0, β1, ϕ, γ, μ, ψ, ω, β0, βreduction1 * β0, βreduction2 * β0)
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection

    sol = memosolver(
        prob, Vern9(; lazy=false); 
        p, u0, callback=cbs, saveat, save_idxs=[ 8 ], 
        abstol=1e-15, maxiters=1e8, verbose=false,
    )
    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return nothing
    end

    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000 .* detection

    for i ∈ eachindex(incidentcases)
        incidence[i] ~ Poisson(incidentcases[i] + 1e-10)
    end
end
