
using DrWatson
@quickactivate :ImmuneBoostingODEs

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple models without having a large 
# number of values in global scope potentially leading to unexpected results

simparms = let 
    γ = 48.7  # generation time 7.5 days
    μ = 0.0087  # Scotland's birth rate = 48000 / 5.5e6
    psis = collect(0:0.1:20)
    R0 = 1.6
    @ntuple γ μ psis R0
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with constant transmission parameter for Fourier transform 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Duration of immunity 6 months 
unforced6mfreqs, unforced6mdensities = simulatedfourierhmdata(0.5, simparms)

let 
    r = findfirst(x -> x >= 1, unforced6mfreqs) # row of `unforced6mdensities` with a period of 1 year
    i = findmax(unforced6mdensities[r, :])
    println("The unforced model has an annual period when ψ = $(collect(0:.1:20)[i[2]])")
end

# Duration of immunity 400 days
unforced400dfreqs, unforced400ddensities = simulatedfourierhmdata((400 / 365), simparms)

# Duration of immunity 2.5 years
unforced25freqs, unforced25densities = simulatedfourierhmdata(2.5, simparms)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with seasonally forced transmission parameter for Fourier transform 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Duration of immunity 6 months 
forced6mfreqs, forced6mdensities, forced6mmins, forced6mmaxs = let 
    @unpack freq_overall, densities = simulatedfourierhmdata(0.5, simparms, 0.1)

    # Calculate values for bifurcation plot 
    @unpack psis = simparms
    mins = Vector{Float64}(undef, length(psis))
    maxs = Vector{Float64}(undef, length(psis))
    for i ∈ eachindex(psis) 
        mins[i] = minimum(results[i].incidence)
        maxs[i] = maximum(results[i].incidence)
    end 
    (freq_overall, densities, mins, maxs,)
end 

# Duration of immunity 400 days
forced400dfreqs, forced400ddensities, forced400dmins, forced400dmaxs = let 
    @unpack freq_overall, densities = simulatedfourierhmdata((400 / 365), simparms, 0.1)

    # Calculate values for bifurcation plot 
    @unpack psis = simparms
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
    @unpack freq_overall, densities = simulatedfourierhmdata(2.5, simparms, 0.1)

    # Calculate values for bifurcation plot 
    @unpack psis = simparms
    mins = Vector{Float64}(undef, length(psis))
    maxs = Vector{Float64}(undef, length(psis))
    for i ∈ eachindex(psis) 
        mins[i] = minimum(results[i].incidence)
        maxs[i] = maximum(results[i].incidence)
    end 
    ( freq_overall, densities, mins, maxs )
end 
