
_tranformmodelparameterproportion(x) = exp(x) / (1 + exp(x))

_inversetranformmodelparameterproportion(x) = log(x / (1 - x))

function _tranformmodelparameters(v; γ, μ, ω)
    @assert length(v) == 8 
    R0 = exp(v[1] + 0.693) 
    β0 = R0 * (γ + μ)
    p = SirnsParameters(
        β0, 
        _tranformmodelparameterproportion(1.5 * v[2]),  # β1 
        _tranformmodelparameterproportion(v[3]) * 2π - π,  # ϕ 
        γ, 
        μ, 
        exp(v[4]),  # ψ
        ω, 
        β0, 
        β0 * _tranformmodelparameterproportion(v[5] + 1.386),  # reducedβ0
        β0 * _tranformmodelparameterproportion(v[6] + 1.386),  # restoredβ0
    )
    pparameter = _tranformmodelparameterproportion(1.5 * v[7])
    detection = _tranformmodelparameterproportion(v[8] - 4.185)
    return ( p, pparameter, detection )
end

@model function fitmodel(
    incidence, prob, cbs, saveat;
    omega,
    rzero_tprior=TDist(2),
    betaone_tprior=TDist(2),
    phi_tprior=TDist(2),
    psi_tprior=TDist(2),
    betareduction1_tprior=TDist(2),
    betareduction2_tprior=TDist(2),
    pparameter_tprior=TDist(2),
    detection_tprior=TDist(2),
)
    rzero_t ~ rzero_tprior
    betaone_t ~ betaone_tprior
    phi_t ~ phi_tprior
    psi_t ~ psi_tprior
    betareduction1_t ~ betareduction1_tprior
    betareduction2_t ~ betareduction2_tprior
    pparameter_t ~ pparameter_tprior
    detection_t ~ detection_tprior

    γ = 48.7  # generation time 7.5 days
    μ = 0.0087  # Scotland's annual birth rate = 48000 / 5.5e6
    ω = omega
    p, pparameter, detection = _tranformmodelparameters(
        ( 
            rzero_t, 
            betaone_t, 
            phi_t, 
            psi_t, 
            betareduction1_t, 
            betareduction2_t, 
            pparameter_t, 
            detection_t 
        ); 
        γ, μ, ω
    )
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection

    sol = solve(
        prob, Vern9(; lazy=false); 
        p, 
        u0, 
        callback=cbs, 
        saveat, 
        save_idxs=[ 8 ], 
        abstol=1e-15, 
        maxiters=5e7,#1e8, 
        verbose=false,
    )
    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return nothing
    end

    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000

    for i ∈ eachindex(incidentcases)
        if isnan(incidentcases[i] * detection * pparameter / (1 - pparameter) + 1e-10)
            Turing.@addlogprob! -Inf
            return nothing
        end
        incidence[i] ~ NegativeBinomial(
            incidentcases[i] * detection * pparameter / (1 - pparameter) + 1e-10,  # > 0
            pparameter
        )
    end
end
