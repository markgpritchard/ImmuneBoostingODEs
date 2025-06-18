
_logistic(x) = 1 / (1 + exp(-x))

@model function fitmodel(
    incidence, prob, cbs, saveat;
    betaone_logitprior=Normal(-2.3, 1.5),  # mean of 0.1 when transformed
    phi_prior=truncated(Normal(0, 1), -π, π),
    gamma_logprior=Normal(log(48.7), 0.4),
    psi_logprior=Normal(0, 2),
    betaprimemultiplier_logitprior=Normal(0, 1.5),
    finalbetaprime_logitprior=Normal(log(0.9 / 0.1), 1.0),
    proportiondetected_logitprior=Normal(log(0.01), 1.0),
    #invrparameter_prior=Exponential(1),
    rparameter_prior=Exponential(1),
    S0_logprior=TDist(2),
    I0_logprior=TDist(2),
    R10_logprior=TDist(2),
    R20_logprior=TDist(2),
    mu=0.0087,
    r0,
    omega,
)
    logitβ1 ~ betaone_logitprior
    ϕ ~ phi_prior
    logγ ~ gamma_logprior
    logψ ~ psi_logprior
    logitbetaprimemultiplier ~ betaprimemultiplier_logitprior
    logitfinalbetaprime ~ finalbetaprime_logitprior
    logitproportiondetected ~ proportiondetected_logitprior
    rparameter ~ rparameter_prior
    logS0 ~ S0_logprior
    logI0 ~ I0_logprior
    logR10 ~ R10_logprior
    logR20 ~ R20_logprior

    # avoid errors caused by extreme values of these parameters
    if isnan(_logistic(logitβ1)) || 
        isnan(_logistic(logitbetaprimemultiplier)) ||
        isnan(_logistic(logitfinalbetaprime)) ||
        isnan(_logistic(logitproportiondetected)) 
        
        Turing.@addlogprob! -Inf
        return nothing
    end

    T = typeof(logitbetaprimemultiplier)
    p = SirnsParameters(
        T(r0 * (exp(logγ) + mu)),  # β0::S
        _logistic(logitβ1),  # β1::T
        ϕ,  # ϕ::T
        exp(logγ),  # γ::T
        mu,  # μ::Float64 
        exp(logψ),  # ψ::T
        omega,  # ω::U
        r0 * (exp(logγ) + mu),  # originalβ0::V
        _logistic(logitbetaprimemultiplier),  # betaprimemultiplier::T
        _logistic(logitfinalbetaprime),  # finalbetaprime::T
        _logistic(logitproportiondetected),  # proportiondetected::T
    )

    _popdenom = 1 + exp(logS0) + exp(logI0) + exp(logR10) + exp(logR20)
    S0 = exp(logS0) / _popdenom
    I0 = exp(logI0) / _popdenom
    R10 = exp(logR10) / _popdenom 
    R20 = exp(logR20) / _popdenom 
    R30 = 1 / _popdenom
    u0 = sirns_u0(S0, I0, R10, R20, R30; p, t0=saveat[1])

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

    cumulativecases = modelcompartments(sol, 1)
    incidentcases = casespertimeblock(cumulativecases) .* 5_450_000

    pparameter = rparameter ./ (rparameter .+ incidentcases .* p.proportiondetected) 

    if rparameter <= 0 || 
        isnan(minimum(pparameter)) || 
        minimum(pparameter) <= 0 || 
        maximum(pparameter) > 1 

        Turing.@addlogprob! -Inf
        return nothing
    end

    incidence ~ arraydist(NegativeBinomial.(rparameter, pparameter))
end
