
using DrWatson
@quickactivate :ImmuneBoostingODEs

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple models without having a large 
# number of values in global scope potentially leading to unexpected results

simparms = let 
    γ = 48.7    # generation time 7.5 days
    μ = 0.0087   # Scotland's birth rate = 48000 / 5.5e6
    psis = collect(0:0.1:20)
    R0 = 1.6
    @ntuple γ μ psis R0
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with constant transmission parameter for Fourier transform 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Duration of immunity 6 months 
unforced6mfreqs, unforced6mdensities = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 0.5
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = ϕ = 0 
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
        )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    fourierhmdata(results)
end 

let 
    # row of `unforced6mdensities` with a period of 1 year
    r = findfirst(x -> x >= 1, unforced6mfreqs) 
    i = findmax(unforced6mdensities[r, :])
    println("The unforced model has an annual period when ψ = $(collect(0:.1:20)[i[2]])")
end

# Duration of immunity 1 year
unforced365dfreqs, unforced365ddensities = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 1
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = ϕ = 0.0  
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
        )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    fourierhmdata(results)
end 

# Duration of immunity 400 days
unforced400dfreqs, unforced400ddensities = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 400 / 365.25
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = ϕ = 0.0  
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
        )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    fourierhmdata(results)
end 

# Duration of immunity 2.5 years
unforced25freqs, unforced25densities = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 2.5
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = ϕ = 0.0  
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
            )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    fourierhmdata(results)
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with seasonally forced transmission parameter for Fourier transform 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Duration of immunity 6 months 
forced6mfreqs, forced6mdensities, forced6mmins, forced6mmaxs = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = .5
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = 0.1
    ϕ = 0.0 
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
        )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    @unpack freq_overall, densities = fourierhmdata(results)

    # Calculate values for bifurcation plot 
    mins = Vector{Float64}(undef, length(psis))
    maxs = Vector{Float64}(undef, length(psis))
    for i ∈ eachindex(psis) 
        mins[i] = minimum(results[i].incidence)
        maxs[i] = maximum(results[i].incidence)
    end 
    ( freq_overall, densities, mins, maxs )
end 

# Duration of immunity 1 year
forced365dfreqs, forced365ddensities, forced365dmins, forced365dmaxs = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 1
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = 0.1
    ϕ = 0.0 
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
        )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    @unpack freq_overall, densities = fourierhmdata(results)

    # Calculate values for bifurcation plot 
    mins = Vector{Float64}(undef, length(psis))
    maxs = Vector{Float64}(undef, length(psis))
    for i ∈ eachindex(psis) 
        mins[i] = minimum(results[i].incidence)
        maxs[i] = maximum(results[i].incidence)
    end 
    ( freq_overall, densities, mins, maxs )
end 

# Duration of immunity 400 days
forced400dfreqs, forced400ddensities, forced400dmins, forced400dmaxs = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 400 / 365.25
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = 0.1
    ϕ = 0.0 
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
        )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    @unpack freq_overall, densities = fourierhmdata(results)

    # Calculate values for bifurcation plot 
    mins = Vector{Float64}(undef, length(psis))
    maxs = Vector{Float64}(undef, length(psis))
    for i ∈ eachindex(psis) 
        mins[i] = minimum(results[i].incidence)
        maxs[i] = maximum(results[i].incidence)
    end 
    ( freq_overall, densities, mins, maxs )
end 

# Duration of immunity 2.5 years
forced25freqs, forced25densities, forced25mins, forced25maxs = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 2.5
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = 0.1
    ϕ = 0.0 
    tspan = ( -1000.0, 20.0 )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(
        undef, length(psis)
    )
    Threads.@threads for i ∈ eachindex(psis)
        ψ = psis[i]
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(
            pl_modelincidence, config, datadir("sims"); 
            prefix="pl_modelincidence"
        )
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    @unpack freq_overall, densities = fourierhmdata(results)

    # Calculate values for bifurcation plot 
    mins = Vector{Float64}(undef, length(psis))
    maxs = Vector{Float64}(undef, length(psis))
    for i ∈ eachindex(psis) 
        mins[i] = minimum(results[i].incidence)
        maxs[i] = maximum(results[i].incidence)
    end 
    ( freq_overall, densities, mins, maxs )
end 
