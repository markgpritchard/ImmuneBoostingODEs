
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version of sirns! with some parameters submitted separately 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

opt_sirns!(du, u, p, t, R0, β1, ψ, fixedbeta1 = nothing) = _opt_sirns!(du, u, p, t, R0, β1, ψ, fixedbeta1)
opt_sirns!(du, u, p, t; R0, β1, ψ, fixedbeta1 = nothing) = opt_sirns!(du, u, p, t, R0, β1, ψ, fixedbeta1)

function _opt_sirns!(du, u, p, t, R0, β1, ψ, fixedbeta1::Nothing)
    βreduction1, βreduction2, βindex = p 
    _opt_sirns!(du, u, p, t, R0, β1, ψ, βreduction1)
end 

function _opt_sirns!(du, u, p, t, R0, β1, ψ, fixedbeta1)
    βreduction1, βreduction2, βindex = p 
    # These parameters are the values used elsewhere in rsvanalysis.j;
    ω = 2
    γ = 48.7    # generation time 7.5 days
    μ = .0087   # Scotland's birth rate = 48000 / 5.5e6
    ϕ = 0
    β0 = R0 * (γ + μ) * [ fixedbeta1, βreduction2 ][round(Int, βindex)]
    newp = @ntuple β0 β1 ϕ ω γ μ ψ
    sirns!(du, u, newp, t)
end 

# run simulation for approx 1000 years to find pre-lockdown stable values to start optimization simulation
function opt_u0(β0, β1, ϕ, γ, μ, ψ, ω, reducetime)
    tspan = ( 1015.35, reducetime )
    p = @ntuple β0 β1 ϕ γ μ ψ ω
    _u0 = sirns_u0(.5, .001; equalrs = true, p, t0 = 1015.35)
    #_u0 = [ [ .5, .001 ]; ones(3) * (1 - .501) / 3; [ cos(2π * .35 - ϕ), sin(2π * .35 - ϕ), 0 ] ]
    prob = ODEProblem(sirns!, _u0, tspan, p)

    sol = solve(prob, Vern9(lazy = false); abstol = 1e-15, reltol = 1e-15, maxiters = 1e7)
    u0 = last(sol)
    u0[8] = 0
    return u0 
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Callbacks for optimization 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# simulate the effect of non-pharmaceutical interventions
opt_changetransmission!(integrator) = integrator.p[3] += 1
# restore β0 
opt_restoretransmission!(integrator) = integrator.p[3] = one(integrator.p[1])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loss function 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

huberloss(a, δ) = a <= δ ? .5 * abs2(a) : δ * (abs(a) - .5 * δ)

function opt_loss(p, prob; casesvector, opt_cbs, saveinds, savetimes, δ = 300, verbose = false)
    saveat = savetimes[saveinds]
    puse = [ p; [ one(p[1]) ] ]
    tmp_prob = remake(prob; p = puse)
    tmp_sol = solve(tmp_prob, Vern9(lazy = false); 
        callback = opt_cbs, saveat, abstol = 1e-12, maxiters = 1e5, reltol = 1e-12)
    if SciMLBase.successful_retcode(tmp_sol.retcode)
        tmp_vector = Array(tmp_sol)[8, :] 
        tmp_casesvector = casesvector[saveinds]
        loss = mean([ 
            huberloss(.0025 * 5.5e6 * (tmp_vector[i+1] - tmp_vector[i]) - case, δ) 
            for (i, case) ∈ enumerate(tmp_casesvector[2:end])
        ])
    else # differential equations solver has not produced sufficiently long simulation
        loss = Inf
    end
    if verbose println(loss) end
    return loss
end

function opt_prob(sirnsfunction, initialvalues, u0, tspan; casesvector, opt_cbs, savetimes, δ = 300)
    saveinds = findall(x -> tspan[1] <= x <= tspan[2], savetimes)
    saveat = savetimes[saveinds]
    p = initialvalues
    prob = ODEProblem(sirnsfunction, u0, tspan, p)
    _loss(p) = opt_loss(p, prob; casesvector, opt_cbs, saveinds, savetimes, δ)
    p0 = initialvalues

    adtype = Optimization.AutoZygote()
    optf = Optimization.OptimizationFunction((x, p) -> _loss(x), adtype)
    optprob = Optimization.OptimizationProblem(optf, p0)

    return optprob
end 

function iterativeopt(R0, β1, ψ, initialvalues; casesvector, increaseday, opt_cbs, reduceday, savetimes) 
    # static values 
    ϕ = -.5π
    γ = 48.7 
    μ = .0087
    immuneduration = .5
    ω = 1 / immuneduration

    β0 = R0 * (γ + μ)
    u0 = opt_u0(β0, β1, ϕ, γ, μ, ψ, ω, reduceday)

    p0 = initialvalues 
    endday = increaseday + .05 
    _sirns1!(du, u, p, t) = opt_sirns!(du, u, p, t; R0, β1, ψ)
    tspan = ( reduceday, endday )
    optprob = opt_prob(_sirns1!, p0, u0, tspan; casesvector, opt_cbs, savetimes, δ = 1000)
    multipliers = Optimization.solve(optprob, ADAM(); maxiters = 500)
    _sirns2!(du, u, p, t) = opt_sirns!(du, u, p, t; R0, β1, ψ, fixedbeta1 = multipliers[1])
    endday += .2 
    p0 = multipliers.u 
    
    while endday < last(savetimes)
        tspan = ( reduceday, endday )
        optprob = opt_prob(_sirns2!, p0, u0, tspan; casesvector, opt_cbs, savetimes)
        multipliers = Optimization.solve(optprob, ADAM(); maxiters = 500)
        endday += .2 
        p0 = multipliers.u 
    end 

    printeffect(p0, R0, β1, ψ)

    return p0
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Display estimated effect of non-pharmaceutical interventions 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function printeffect(multipliers, R0, β1, ψ)
    reduction = round(1 - multipliers[1]; sigdigits = 3) 
    restore = round(multipliers[2]; sigdigits = 3)
    intro = "For R0 = $R0, β1 = $β1, ψ = $ψ"
    initialeffect = "estimated reduction in transmission from non-pharmaceutical interventions is"
    latereffect = "and after interventions, transmission returned to"
    println("$intro, $initialeffect $reduction, $latereffect $restore of the original")
end 
