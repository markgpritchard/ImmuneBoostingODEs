
_logistic(x) = exp(x) / (1 + exp(x))
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
    rzero_prior=Gamma(1, 2),
    betaone_logitprior=Normal(0, 2.0),
    phi_prior=Uniform(-π, π),
    psi_logprior=Normal(0, 2.0),
    omega_logprior=Normal(0, 2.0),
    betaprimemultiplier_logitprior=Normal(0, 2.0),
    finalbetaprime_logitprior=Normal(0, 2.0),
    proportiondetected_logitprior=Normal(0, 2.0),
    #pparameter_tranformedlogitprior=TDist(2),
    #pparameter_tranformedlogitprior=truncated(TDist(2), -2.1, Inf),
    #pparameter_logitprior=truncated(TDist(2), -1.0, Inf),
    #pparameter_logitprior=truncated(Normal(0, 1), -1.0, Inf),
    #inverserprior=Exponential(1),
    rparameter_prior=Exponential(1),
    gamma=48.7,
    mu=0.0087,
)
    r0 ~ rzero_prior
    logitβ1 ~ betaone_logitprior
    ϕ ~ phi_prior
    logψ ~ psi_logprior
    logω ~ omega_logprior
    logitbetaprimemultiplier ~ betaprimemultiplier_logitprior
    logitfinalbetaprime ~ finalbetaprime_logitprior
    logitproportiondetected ~ proportiondetected_logitprior
    #tranformedlogitpparameter ~ pparameter_tranformedlogitprior
    #logitpparameter ~ pparameter_logitprior
    rparameter ~ rparameter_prior

    if isnan(r0) || rparameter <= 0
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
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection

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
