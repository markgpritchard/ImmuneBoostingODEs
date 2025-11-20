
function fitsirns_u0(S0, I0, R1, R2; p, t0=0)
    # accepts values summing to greater than 1 and standardizes
    tot = sum([S0, I0, R1, R2, 0.1])
    return sirns_u0(S0 / tot, I0 / tot, R1 / tot, R2 / tot, 0.1 / tot; p, t0)
end

function sirnsloss(
    parms;  
    incidence, 
    prob,
    callback,
    saveat,
    gamma=48.7,  # generation time 7.5 days
    mu=0.0087,  # Scotland's birth rate = 48000 / 5.5e6
    psi=0,
    population=5_450_000,
    abstol=1e-15, 
    maxiters=1e8,    
)
    p = SirnsParameters( ; 
        β0=exp(parms[1]), 
        β1=exp(parms[2]), 
        ϕ=exp(parms[3]), 
        γ=gamma, 
        μ=mu, 
        ψ=psi, 
        ω=exp(parms[4]), 
        βreductionfactor=exp(parms[5]), 
        βreduction=1.0,
    )
    detection = exp(parms[6])
    minsigma2 = exp(parms[7])
    u0 = fitsirns_u0(exp.(parms[8:11])...; p, t0=2005.737)  # 1 years before data collection
    

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
    if sol.retcode != :Success
        return Inf
    end

    cumulativecases = modelcompartments(sol, 1)
    weeklyincidentcases = casespertimeblock(cumulativecases) .* population

    np = weeklyincidentcases .* detection
    npminus = np .* (1 - detection) .+ minsigma2
    
    # Normal approximation of Binomial to avoid forcing integer values 
    #incidence ~ arraydist(Normal.(np, NaNMath.sqrt.(npminus)))=#

    return sum(-log.(pdf.(Normal.(np, NaNMath.sqrt.(npminus)), incidence)))
end

function optimizesirns(
    incidence, 
    prob, 
    parms; # which is a vector containing the log of the following parameters in order:
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
    callback,
    saveat,
    gamma=48.7,  # generation time 7.5 days
    mu=0.0087,  # Scotland's birth rate = 48000 / 5.5e6
    psi=0,
    population=5_450_000,
    abstol=1e-15, 
    maxiters=1e8, 
    lb, 
    ub,
    optimizationsolvermaxiters=1e5,
    adtype=Optimization.AutoZygote(),
    nt=10,
    rt=0.975,
    r_expand=2.0,
    verbosity=3,
)
    p = SirnsParameters( ; 
        β0=exp(parms[1]), 
        β1=exp(parms[2]), 
        ϕ=exp(parms[3]), 
        γ=gamma, 
        μ=mu, 
        ψ=psi, 
        ω=exp(parms[4]), 
        βreductionfactor=exp(parms[5]), 
        βreduction=1.0,
    )
    u0 = fitsirns_u0(exp.(parms[6:9])...; p, t0=2005.737)  # 1 years before data collection
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
        @error "$(sol.retcode) with parms=$parms, u0=$u0"
        return nothing
    end
    
    optf = Optimization.OptimizationFunction(
        (x, p) -> sirnsloss(
            x; 
            incidence, 
            prob,
            callback,
            saveat,
            gamma,
            mu,
            psi,
            population,
            abstol, 
            maxiters,  
        ), 
        adtype
    )
    optprob = Optimization.OptimizationProblem(optf, parms; lb, ub)
    result_ode = Optimization.solve(
        optprob, Optim.SAMIN(; nt, rt, r_expand, verbosity); 
        maxiters=optimizationsolvermaxiters
    )

    return result_ode
end

@model function fitmodel(
    incidence, prob;
    callback,
    saveat,
    betazeroprior=Exponential(100), 
    betaoneprior=Uniform(0, 0.9),
    phiprior=Uniform(-π, π),
    gamma=48.7,  # generation time 7.5 days
    mu=0.0087,  # Scotland's birth rate = 48000 / 5.5e6
    psi=0,
    omegaprior=LogNormal(0, 1),
    betareductionfactorprior=LogNormal(0, 1),
    S0prior=LogNormal(log(0.5), 1),
    I0prior=LogNormal(log(0.1), 1),
    R10prior=LogNormal(log(0.1), 1),
    R20prior=LogNormal(log(0.1), 1),
    detectionprior=Beta(1, 98),
    population=5_450_000,
    abstol=1e-15, 
    maxiters=1e8, 
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
    minsigma2 ~ Beta(1, 2)
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
