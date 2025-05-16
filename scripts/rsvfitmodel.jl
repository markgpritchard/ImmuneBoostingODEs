
_logistic(x) = 1 / (1 + exp(-x))
#=
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
=#
@model function fitmodel(
    incidence, prob, cbs, saveat;
    #rzero_prior=Gamma(2, 1),
    #rzero_prior=Gamma(1.5, 1),
    rzero_logprior=Normal(0.7, 1.2),
    betaone_logitprior=Normal(-2.3, 1.5),  # mean of 0.1
    phi_prior=truncated(Normal(0, 1), -π, π),
    psi_logprior=Normal(0, 2),
    #psi_logprior=TDist(2),
    omega_logprior=Normal(0, 0.7),
    betaprimemultiplier_logitprior=Normal(0, 1.5),
    finalbetaprime_logitprior=Normal(log(0.9 / 0.1), 1.0),
    proportiondetected_logitprior=Normal(log(0.01), 1.0),
    #pparameter_tranformedlogitprior=TDist(2),
    #pparameter_tranformedlogitprior=truncated(TDist(2), -2.1, Inf),
    #pparameter_logitprior=truncated(TDist(2), -1.0, Inf),
    #pparameter_logitprior=truncated(Normal(0, 1), -1.0, Inf),
    #inverserprior=Exponential(1),
    invrparameter_prior=Exponential(1),
    S0max_logitprior=TDist(2),
    I0_transformedlogitprior=TDist(2),
    gamma=48.7,
    mu=0.0087,
)
    logr0 ~ rzero_logprior
    logitβ1 ~ betaone_logitprior
    ϕ ~ phi_prior
    logψ ~ psi_logprior
    logω ~ omega_logprior
    logitbetaprimemultiplier ~ betaprimemultiplier_logitprior
    logitfinalbetaprime ~ finalbetaprime_logitprior
    logitproportiondetected ~ proportiondetected_logitprior
    #tranformedlogitpparameter ~ pparameter_tranformedlogitprior
    #logitpparameter ~ pparameter_logitprior
    invrparameter ~ invrparameter_prior
    logitS0max ~ S0max_logitprior
    transformedlogitI0 ~ I0_transformedlogitprior

    r0 = exp(logr0)

    if isnan(r0) || 1 / invrparameter <= 0
        Turing.@addlogprob! -Inf
        return nothing
    end

    p = SirnsParameters(
        r0 * (gamma + mu),  # β0::T
        _logistic(logitβ1),  # β1::T
        ϕ,  # ϕ::T
        gamma,  # γ::Float64
        mu,  # μ::Float64 
        exp(logψ),  # ψ::T
        exp(logω),  # ω::T
        r0 * (gamma + mu),  # originalβ0::T
        _logistic(logitbetaprimemultiplier),  # betaprimemultiplier::T
        _logistic(logitfinalbetaprime),  # finalbetaprime::T
        _logistic(logitproportiondetected),  # proportiondetected::T
    )
#=
    pparameter = _logistic(logitpparameter)
    if pparameter == 0 || pparameter == 1
        Turing.@addlogprob! -Inf
        return nothing
    end
 =#
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
        #=
        if isnan(
            incidentcases[i] * p.proportiondetected * pparameter / (1 - pparameter) + 1e-10
            )
            Turing.@addlogprob! -Inf
            return nothing
        end
        incidence[i] ~ NegativeBinomial(
            incidentcases[i] * p.proportiondetected * pparameter / (1 - pparameter) + 1e-10,  # > 0
            pparameter
        )=#
        #=incidence[i] ~ NegativeBinomial(
            rparameter + 1e-10,
            (rparameter + 1e-10) / ((rparameter + 1e-10) + incidentcases[i] * p.proportiondetected)
        )=#
        incidence[i] ~ NegativeBinomial(
            rparameter,
            rparameter / (rparameter + incidentcases[i] * p.proportiondetected)
        )
    end
end
