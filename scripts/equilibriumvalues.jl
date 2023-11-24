
# Find and plot equilibrium values for model 

using DrWatson
@quickactivate :ImmuneBoostingODEs

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple simulations, and across running 
# and plotting simulations, without having a large number of values in global scope 
# potentially leading to unexpected or inconsistent results

equilparms = let 
    γ       = 48.7    # generation time 7.5 days
    μ       = .0087   # Scotland's birth rate = 48_000 / 5.5e6
    psis    = collect(0:.1:20)
    ω       = .913    # mean duration of immunity without boosting 400 days 
    durations = [ collect(.01:.01:2.5); collect(2.5:.1:80) ] # Used for plotting 
    # effect of changing mean duration of immunity. Set to provide finer resolution 
    # for first 2.5 years
    R0s     = collect(.0125:.025:15)
    @ntuple γ μ psis ω durations R0s
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Equilibrium proportion in each compartment for different parameter values 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Values of I^* for different values of R0 and ϕ, which is plotted in Figure 1

let 
    @unpack γ, μ, psis, ω, R0s = equilparms
    betas = R0s .* (γ + μ)
    
    global equilSs = zeros(length(betas), length(psis))
    global equilIs = zeros(length(betas), length(psis))
    global equilRs = zeros(length(betas), length(psis))
    for (i, β) ∈ enumerate(betas), (j, ϕ) ∈ enumerate(psis)
        I = equili(SirnsParameters(β, γ, μ, ϕ, ω))
        if 0 < I <= 1   # endemic equilibrium
            equilSs[i, j] = equils(SirnsParameters(β, γ, μ, ϕ, ω))
            equilIs[i, j] = I
            equilRs[i, j] = equilr(SirnsParameters(β, γ, μ, ϕ, ω), I) 
        else            # disease-free equilibium
            equilSs[i, j] = 1.
            equilIs[i, j] = .0
            equilRs[i, j] = .0
        end 
    end 
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stability of endemic equilibria for different values of ω 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For Figure 2 we want the stability for a range of values of immune duration at 
# different R0 values

critvector = let 
    @unpack γ, μ, psis, durations = equilparms
    R0s = [ 1.6, 5, 10, 15 ]
    betas = R0s .* (γ + μ)
    omegas = 1 ./ durations
    [ [ realmaxequileigen(SirnsParameters(β, γ, μ, ψ, ω); warntol = Inf) 
            for ψ ∈ psis, ω ∈ omegas ]
        for β ∈ betas ]
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations to show stability and limit cycles
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Examples with R0 = 1.6 are for Figure 1, and examples with other R0 values for 
# supplementary figures.
 
# find the relevant critical value of ϕ 
let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    β = R0 * (γ + μ)
    println("ϕ = $(findpsi(β, γ, μ, ω; warntol = Inf))")
end     # ϕ = 1.18

# stable simulation
stablesim = let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    ϕ = .5
    β = R0 * (γ + μ)
    p = SirnsParameters(β, γ, μ, ϕ, ω)
    u0 = sirns_u0(.8, .03; p, equalrs = true)
    tspan = ( 0., 35. ) # it has to run a long time to avoid looking like a limit cycle
    sol = run_sirns(u0, p, tspan)
    modelcompartments(sol, p)
end

# unstable simulations
let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    ϕ = 1.5
    β = R0 * (γ + μ)
    p = SirnsParameters(β, γ, μ, ϕ, ω)
    
    # find the limit cycle 
    let
        u0 = sirns_u0(.7, .1; p, equalrs = true)
        tspan = ( -1000., 10. )
        sol = run_sirns(u0, p, tspan)
        global limitcycle12 = modelcompartments(sol, p)
    end

    # approach the cycle from inside
    let
        u0 = sirns_u0(.64, .008; p, equalrs = true)
        tspan = ( 0., 3.8 )
        sol = run_sirns(u0, p, tspan)
        global unstablein12 = modelcompartments(sol, p)
    end

    # approach the cycle from outside
    let
        u0 = sirns_u0(.8, .03; p, equalrs = true)
        tspan = ( 0., 3.1 )
        sol = run_sirns(u0, p, tspan)
        global unstableout12 = modelcompartments(sol, p)
    end 
end 

let
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    ϕ = 5
    β = R0 * (γ + μ)
    p = SirnsParameters(β, γ, μ, ϕ, ω)
    
    # find the limit cycle 
    let
        u0 = sirns_u0(.7, .1; p, equalrs = true)
        tspan = ( -1000., 10. )
        sol = run_sirns(u0, p, tspan)
        global limitcycle5 = modelcompartments(sol, p)
    end

    # approach the cycle from inside
    let
        u0 = sirns_u0(.75, .01; p, equalrs = true)
        tspan = ( 0., 3.2 )
        sol = run_sirns(u0, p, tspan)
        global unstablein5 = modelcompartments(sol, p)
    end

    # approach the cycle from outside
    let
        u0 = sirns_u0(.8, .03; p, equalrs = true)
        tspan = ( 0., 3.5 )
        sol = run_sirns(u0, p, tspan)
        global unstableout5 = modelcompartments(sol, p)
    end 
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

let 
    @unpack γ, μ, ω = equilparms
    R0 = 1.6
    β = R0 * (γ + μ)
    tspan = ( 0., 10. )

    let 
        ϕ = .5 
        p = SirnsParameters(β, γ, μ, ϕ, ω)
        u0 = sirns_u0(.7, .01; p, equalrs = true)
        sol = run_sirns(u0, p, tspan)
        global basic5 = modelcompartments(sol, p)
    end 

    let 
        ϕ = 1.5 
        p = SirnsParameters(β, γ, μ, ϕ, ω)
        u0 = sirns_u0(.7, .01; p, equalrs = true)
        sol = run_sirns(u0, p, tspan)
        global basic1_5 = modelcompartments(sol, p)
    end 

    let 
        ϕ = 5 
        p = SirnsParameters(β, γ, μ, ϕ, ω)
        u0 = sirns_u0(.7, .01; p, equalrs = true)
        sol = run_sirns(u0, p, tspan)
        global basic5_ = modelcompartments(sol, p)
    end 
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Values for bifurcation plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bifurcationI_1_5 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 1.5
    config = @dict R0 γ μ psis ω
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end

bifurcationI_5 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 5
    config = @dict R0 γ μ psis ω maxiters = 1e6
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end

bifurcationI_10 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 10
    config = @dict R0 γ μ psis ω maxiters = 1e6
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end

bifurcationI_15 = let
    @unpack γ, μ, psis, ω = equilparms
    R0 = 15
    config = @dict R0 γ μ psis ω maxiters = 1e6
    d = produce_or_load(pl_bifurcationlimits, config, datadir("sims"); 
        prefix = "pl_bifurcationlimits")
    dict2ntuple(d[1])
end
