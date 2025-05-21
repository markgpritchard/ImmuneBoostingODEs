
_logistic(x) = 1 / (1 + exp(-x))

@model function fitmodel(
    incidence, prob, cbs, saveat;
    betaone_logitprior=Normal(-2.3, 1.5),  # mean of 0.1 when transformed
    phi_prior=truncated(Normal(0, 1), -π, π),
    psi_logprior=Normal(0, 2),
    betaprimemultiplier_logitprior=Normal(0, 1.5),
    finalbetaprime_logitprior=Normal(log(0.9 / 0.1), 1.0),
    proportiondetected_logitprior=Normal(log(0.01), 1.0),
    invrparameter_prior=Exponential(1),
    S0max_logitprior=TDist(2),
    I0_transformedlogitprior=TDist(2),
    gamma=48.7,
    mu=0.0087,
    r0,
    omega,
)
    logitβ1 ~ betaone_logitprior
    ϕ ~ phi_prior
    logψ ~ psi_logprior
    logitbetaprimemultiplier ~ betaprimemultiplier_logitprior
    logitfinalbetaprime ~ finalbetaprime_logitprior
    logitproportiondetected ~ proportiondetected_logitprior
    invrparameter ~ invrparameter_prior
    logitS0max ~ S0max_logitprior
    transformedlogitI0 ~ I0_transformedlogitprior

    if 1 / invrparameter <= 0
        Turing.@addlogprob! -Inf
        return nothing
    end

    T = typeof(logitbetaprimemultiplier)
    p = SirnsParameters(
        T(r0 * (gamma + mu)),  # β0::S
        _logistic(logitβ1),  # β1::T
        ϕ,  # ϕ::T
        gamma,  # γ::Float64
        mu,  # μ::Float64 
        exp(logψ),  # ψ::T
        omega,  # ω::U
        r0 * (gamma + mu),  # originalβ0::V
        _logistic(logitbetaprimemultiplier),  # betaprimemultiplier::T
        _logistic(logitfinalbetaprime),  # finalbetaprime::T
        _logistic(logitproportiondetected),  # proportiondetected::T
    )

    I0 = _logistic(transformedlogitI0 - 6)
    S0 = min(_logistic(logitS0max), 1 - I0)
    u0 = sirns_u0(S0, I0; p, equalrs=true, t0=1996.737)  # 10 years before data collection

    sol = solve(
        prob, Vern9(; lazy=false); 
        p, 
        u0, 
        callback=cbs, 
        saveat, 
        save_idxs=[ 8 ], 
        abstol=1e-15, 
        maxiters=1e8, 
        verbose=false,
    )
    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return nothing
    end

    rparameter = 1 / invrparameter
    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000

    for i ∈ eachindex(incidentcases)
        incidence[i] ~ NegativeBinomial(
            rparameter,
            rparameter / (rparameter + incidentcases[i] * p.proportiondetected)
        )
    end
end
