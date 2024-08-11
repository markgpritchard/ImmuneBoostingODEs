
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fourier analysis of frequency of oscillations 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function fouriervalues(modeloutput)
    @unpack gt, incidence = modeloutput
    return fouriervalues(incidence, gt)
end 

function fouriervalues(values::Vector, t::Vector)
    dt = t[2] - t[1]
    n = length(t) 
    return fouriervalues(values, dt, n)
end 

function fouriervalues(values, dt::Real, n::Int)
    stdvalues = (values .- mean(values)) ./ std(values)
    valueshat = fft(stdvalues)
    spectraldensity = @. abs(valueshat) / n
    freq = 1 / (dt * n) * (1:n)
    inds = Int.(1:floor(n / 2))
    return (freq = freq[inds], spectraldensity = spectraldensity[inds])
end 

"""
    fourierhmdata(modeloutputs)

Performs fast Fourier transform and returns a `NamedTuple` of frequencies and densities.

`modeloutputs` is a `NamedTuple` containing `gt`, times points from a simulation, 
    and `incidence`, a vector of incidence from a simulation.
"""
function fourierhmdata(modeloutputs) 
    @unpack freq, spectraldensity = fouriervalues(modeloutputs[1])
    freq_overall = deepcopy(freq)
    densities = zeros(length(spectraldensity), length(modeloutputs))
    fourierhmdata!(densities, modeloutputs, freq_overall)
    return @ntuple freq_overall densities
end 

function fourierhmdata!(densities, modeloutputs, freq_overall)
    for (i, m) ∈ enumerate(modeloutputs) 
        fourierhmdata!(densities, m, freq_overall, i) 
    end 
end 

function fourierhmdata!(densities, m, freq_overall, i)
    @unpack freq, spectraldensity = fouriervalues(m)
    @assert freq == freq_overall 
    densities[:, i] = spectraldensity 
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter fitting 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

memosolver(prob; kwargs...) = @memoize solve(prob; kwargs...)
memosolver(prob, alg; kwargs...) = @memoize solve(prob, alg; kwargs...)

function loadrsvdata(omega; ids=1:4, maxrounds=12,)
    df = DataFrame(
        :iteration => Int[ ],
        :chain => Int[ ],
        :β0 => Float64[ ],
        :β1 => Float64[ ],
        :ϕ => Float64[ ],
        :ψ => Float64[ ],
        :βreduction1 => Float64[ ],
        :βreduction2 => Float64[ ],
        :log_density => Float64[ ],
    )
    for i ∈ ids
        _loaded = false 
        j = maxrounds
        while !_loaded && j >= 1
            filename = "rsvparameters_omega_$(omega)_id_$(i)_nrounds_$(j).jld2"
            if isfile(datadir("sims", filename))
                _df = DataFrame(load(datadir("sims", filename))["chain"])
                _df.chain = [ i for _ ∈ axes(_df, 1) ]
                df = vcat(df, _df)
                _loaded = true
            end
            j += -1
        end
    end
    return df
end

function fittedsimulationsetup(saveat)
    p = SirnsParameters(50.0, 0.1, 0.0, 48.7, 0.0087, 0.0, 2.0, 50.0, 50.0, 50.0) 
    u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)
    tspan = ( 1996.737, last(saveat) )
    ODEProblem(sirns!, u0, tspan, p)
end

function runfittedsimulations(df, omega, saveat, cbs; detection=0.02)
    modelmat = Matrix{Float64}(undef, length(saveat) - 1, size(df, 1))
    prob = fittedsimulationsetup(saveat)
    for i ∈ axes(df, 1)
        p = SirnsParameters(
            df.β0[i], 
            df.β1[i], 
            df.ϕ[i], 
            48.7, 
            0.0087, 
            df.ψ[i], 
            omega, 
            df.β0[i], 
            df.βreduction1[i] * df.β0[i], 
            df.βreduction2[i] * df.β0[i]
        )
        u0 = sirns_u0(0.01, 2e-5; p, equalrs=true, t0=1996.737)
        sol = memosolver(
            prob, Vern9(; lazy=false); 
            p, u0, callback=cbs, saveat, abstol=1e-15, maxiters=1e8, 
        )
        cumulativecases = modelcompartments(sol, :cc)
        modelmat[:, i] = casespertimeblock(cumulativecases) .* 5_450_000 .* detection
    end
    return modelmat
end

function fittedsimulationquantiles(
    modelmat::AbstractMatrix, quantiles::AbstractVector=[ 0.05, 0.5, 0.95 ]
)
    return [ quantile(modelmat[i, :], quantiles) for i ∈ axes(modelmat, 1) ]
end

function fittedsimulationquantiles(
    df::DataFrame, omega::Number, saveat::AbstractVector, cbs, 
    quantiles::AbstractVector=[ 0.05, 0.5, 0.95 ]
)
    modelmat = runfittedsimulations(df, omega, saveat, cbs)
    return fittedsimulationquantiles(modelmat, quantiles)
end
