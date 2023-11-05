
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fourier analysis of frequency of oscillations 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
