
struct FailedU0 end 
const failedu0 = FailedU0()

_validu0(::FailedU0) = false 
_validu0(::AbstractVector) = true

const DEFAULTPRIORS = [
    truncated(Exponential(100), 0.5 * 48.7, 25 * 48.7),  # β0 
    Beta(1, 1),  # β1 
    Uniform(-π, π),  # ϕ 
    LogNormal(0, 1),  # ω 
    LogNormal(0, 1),  # βreductionfactor 
    Beta(1, 98),  # detection 
    Beta(1, 2),  # minsigma2 
    LogNormal(log(0.5), 1),  # S0 
    LogNormal(log(0.1), 1),  # I0 
    LogNormal(log(0.1), 1),  # R1 
    LogNormal(log(0.1), 1),  # R2 
]

function fitsirns_u0(S0, I0, R1, R2; p, t0=0)
    #println("S0=$S0, I0=$I0, R1=$R1, R2=$R2, p=$p")
    # accepts values summing to greater than 1 and standardizes
    tot = sum([S0, I0, R1, R2, 0.1])
    if isnan(tot)
        @warn "NaN values passed to u0, S0=$S0, I0=$I0, R1=$R1, R2=$R2, p=$p"
        return failedu0
    else
        return sirns_u0(S0 / tot, I0 / tot, R1 / tot, R2 / tot, 0.1 / tot; p, t0)
    end
end

function initialparams_map(model, nchains; kwargs...)
    return _initialparams(maximum_a_posteriori, model, nchains; kwargs...)
end

function initialparams_mle(model, nchains; kwargs...)
    return _initialparams(maximum_likelihood, model, nchains; kwargs...)
end

function _initialparams(f, model, nchains; kwargs...)
    initial_params = Vector{InitFromParams}(undef, nchains)
    _initialparams!(f, initial_params, model, nchains; kwargs...)
    return initial_params
end

function _initialparams!(
    f, initial_params, model, nchains; 
    displayparams=true, 
    lb=quantile.(DEFAULTPRIORS, 0.005),
    ub=quantile.(DEFAULTPRIORS, 0.995),
    maxiters=1000, 
    kwargs...
)
    Threads.@threads for k in 1:nchains 
        paramvalues = f(model, LBFGS(; linesearch=Static()); lb, ub, maxiters, kwargs...)
        initial_params[k] = InitFromParams(paramvalues.params)
        if displayparams 
            rc = paramvalues.optim_result.retcode
            @info "initial_params[$k] ($rc) = $(initial_params[k])"
        end
    end

    return nothing
end

@model function fitmodel(
    incidence, prob;
    callback,
    saveat,
    betazeroprior=DEFAULTPRIORS[1], 
    betaoneprior=DEFAULTPRIORS[2], 
    phiprior=DEFAULTPRIORS[3], 
    omegaprior=DEFAULTPRIORS[4], 
    betareductionfactorprior=DEFAULTPRIORS[5], 
    detectionprior=DEFAULTPRIORS[6], 
    minsigma2prior=DEFAULTPRIORS[7], 
    S0prior=DEFAULTPRIORS[8], 
    I0prior=DEFAULTPRIORS[9], 
    R10prior=DEFAULTPRIORS[10], 
    R20prior=DEFAULTPRIORS[11], 
    gamma=48.7,  # generation time 7.5 days
    mu=0.0087,  # Scotland's birth rate = 48000 / 5.5e6
    psi=0,
    population=5_450_000,
    abstol=1e-15, 
    maxiters=5e5,#1e8 
)
    β0 ~ betazeroprior
    β1 ~ betaoneprior
    ϕ ~ phiprior
    γ = gamma  
    μ = mu  # Scotland's birth rate = 48000 / 5.5e6
    ψ = psi
    ω ~ omegaprior
    βreductionfactor ~ betareductionfactorprior
    detection ~ detectionprior
    minsigma2 ~ minsigma2prior
    S0 ~ S0prior
    I0 ~ I0prior
    R1 ~ R10prior
    R2 ~ R20prior

    p = SirnsParameters( ; 
        β0, 
        β1, 
        ϕ, 
        γ, 
        μ, 
        ψ, 
        ω, 
        βreductionfactor, 
        βreduction=1.0,
    )
    #u0 = sirns_u0(min(mu / gamma, 1 - 2e-5), 2e-5; p, equalrs=true, t0=1996.737)  # 10 years before data collection
    u0 = fitsirns_u0(S0, I0, R1, R2; p, t0=2005.737)  # 1 years before data collection

    if !_validu0(u0)
        @addlogprob! -Inf
        return nothing
    end

    #sol = memosolver(
    sol = solve(
        prob, Vern9(; lazy=false); 
        p, 
        u0, 
        callback, 
        saveat, 
        save_idxs=[8], 
        abstol, 
        maxiters, 
        verbose=false,
    )
    if !successful_retcode(sol)
        @addlogprob! -Inf
        return nothing
    end

    cumulativecases = modelcompartments(sol, 1)
    weeklyincidentcases = casespertimeblock(cumulativecases) .* population

    # Normal approximation of Binomial to avoid forcing integer values 
    np = weeklyincidentcases .* detection
    npminus = np .* (1 - detection) .+ minsigma2
    
    # Normal approximation of Binomial to avoid forcing integer values 
    incidence ~ arraydist(Normal.(np, NaNMath.sqrt.(npminus)))
    return nothing
end
