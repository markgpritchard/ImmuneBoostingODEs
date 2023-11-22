
using DrWatson
@quickactivate :ImmuneBoostingODEs
using CairoMakie
CairoMakie.activate!() # allows figures to be saved as vector files

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters used frequently in this script 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Allows the same parameters to be used in multiple models without having a large 
# number of values in global scope potentially leading to unexpected results

simparms = let 
    γ = 48.7    # generation time 7.5 days
    μ = .0087   # Scotland's birth rate = 48000 / 5.5e6
    psis = collect(0:.1:20)
    R0 = 1.6
    @ntuple γ μ psis R0
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with constant transmission parameter for Fourier transform 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Duration of immunity 6 months 
unforced6mfreqs, unforced6mdensities = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = .5
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = ϕ = 0 
    tspan = ( -1000., 20. )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(undef, length(psis))
    for (i, ψ) ∈ enumerate(psis)
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(pl_modelincidence, config, datadir("sims"); 
            prefix = "pl_modelincidence")
        @unpack gt, incidence = pl_result[1]
        results[i] = @ntuple gt incidence
    end
    fourierhmdata(results)
end 

let 
    r = findfirst(x -> x >= 1, unforced6mfreqs) # row of `unforced6mdensities` with a period of 1 year
    i = findmax(unforced6mdensities[r, :])
    println("The unforced model has an annual period when ψ = $(collect(0:.1:20)[i[2]])")
end

# Duration of immunity 400 days
unforced400dfreqs, unforced400ddensities = let 
    @unpack γ, μ, psis, R0 = simparms
    immuneduration = 400 / 365.25
    β0 = R0 * (γ + μ)
    ω = 1 / immuneduration 
    β1 = ϕ = 0  
    tspan = ( -1000., 20. )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(undef, length(psis))
    for (i, ψ) ∈ enumerate(psis)
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(pl_modelincidence, config, datadir("sims"); 
            #filename = hash, 
            prefix = "pl_modelincidence")
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
    β1 = ϕ = 0  
    tspan = ( -1000., 20. )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(undef, length(psis))
    for (i, ψ) ∈ enumerate(psis)
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(pl_modelincidence, config, datadir("sims"); 
            prefix = "pl_modelincidence")
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
    β1 = .1
    ϕ = 0 
    tspan = ( -1000., 20. )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(undef, length(psis))
    for (i, ψ) ∈ enumerate(psis)
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(pl_modelincidence, config, datadir("sims"); 
            prefix = "pl_modelincidence")
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
    β1 = .1
    ϕ = 0 
    tspan = ( -1000., 20. )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(undef, length(psis))
    for (i, ψ) ∈ enumerate(psis)
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(pl_modelincidence, config, datadir("sims"); 
            prefix = "pl_modelincidence")
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
    β1 = .1
    ϕ = 0 
    tspan = ( -1000., 20. )
    kw = @ntuple tspan
    results = Vector{NamedTuple{(:gt, :incidence), Tuple{Vector{Float64}, Vector{Float64}}}}(undef, length(psis))
    for (i, ψ) ∈ enumerate(psis)
        config = @dict β0 β1 ϕ γ μ ψ ω kw
        pl_result = produce_or_load(pl_modelincidence, config, datadir("sims"); 
            prefix = "pl_modelincidence")
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots of simulations 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fourierplot = Figure(; size = ( 800, 300 ))
fourierplot_gb = GridLayout(fourierplot[1, 1])

Label(fourierplot_gb[0, 1], "Immune duration 6 months"; 
    fontsize = 11.84, halign = :left, tellwidth = false)
Label(fourierplot_gb[0, 2], "1.1 years"; fontsize = 11.84, halign = :left, tellwidth = false)
Label(fourierplot_gb[0, 3], "2.5 years"; fontsize = 11.84, halign = :left, tellwidth = false)

let 
    @unpack psis = simparms
    freqvector = [ unforced6mfreqs, unforced400dfreqs, unforced25freqs ]
    densityvector = [ unforced6mdensities, unforced400ddensities, unforced25densities ]

    axs = [ Axis(fourierplot_gb[1, i]; yscale = log) for i ∈ 1:3 ]
    i = 1
    for (ax, freqs, densities) ∈ zip(axs, freqvector, densityvector)
        inds = findall(x -> .1 <= x <= 365 / 21, freqs) 
        hm = heatmap!(ax, psis, freqs[inds], log10.(densities[inds, :]'); 
            colormap = :CMRmap, colorrange = ( -6, 0 ))
        ax.yticks = ([.2, .5, 1, 2, 4, 365 / 35], [ "5 y", "2 y", "1 y", "6 m", "3 m", "5 w"])
        formataxis!(ax; hidex = true, hidey = i != 1)
        if i == 1 
            co = Colorbar(fourierplot_gb[1:2, 4], hm) 
            formataxis!(co)
        end
        i += 1
    end
    Label(fourierplot_gb[1:2, 0], "Period"; 
        fontsize = 11.84, rotation = π/2, tellheight = false)
    Label(fourierplot_gb[1:2, 5], "log₁₀ spectral density"; 
        fontsize = 11.84, rotation = -π/2, tellheight = false)
end

let 
    @unpack psis = simparms
    freqvector = [ forced6mfreqs, forced400dfreqs, forced25freqs ]
    densityvector = [ forced6mdensities, forced400ddensities, forced25densities ]

    axs = [ Axis(fourierplot_gb[2, i]; yscale = log) for i ∈ 1:3 ]
    i = 1
    for (ax, freqs, densities) ∈ zip(axs, freqvector, densityvector)
        inds = findall(x -> .1 <= x <= 365 / 21, freqs) 
        hm = heatmap!(ax, psis, freqs[inds], log10.(densities[inds, :]'); 
            colormap = :CMRmap, colorrange = ( -6, 0 ))
        ax.yticks = ([.2, .5, 1, 2, 4, 365 / 35], [ "5 y", "2 y", "1 y", "6 m", "3 m", "5 w"])
        formataxis!(ax; hidey = i != 1)
        i += 1
    end
    Label(fourierplot_gb[3, 1:3], "Boosting coefficient, ψ"; 
        fontsize = 11.84, tellwidth = false)
end

labelplots!([ "A", "B"], fourierplot_gb; rows = [ 0, 2])

for r ∈ [ 1, 3 ] rowgap!(fourierplot_gb, r, 5) end
for c ∈ [ 1, 5 ] colgap!(fourierplot_gb, c, 5) end
fourierplot

safesave(plotsdir("fourierplot.pdf"), fourierplot)
