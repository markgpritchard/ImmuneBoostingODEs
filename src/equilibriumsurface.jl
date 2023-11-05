
function surfacegrid!(ax, xs, ys, zs; nsamples = 20, linewidth = .25, rasterize = 20, kwargs...)
    xind = surfacegrid_samples(xs, nsamples)
    yind = surfacegrid_samples(ys, nsamples)

    # The plots 
    surface!(ax, xs, ys, zs; rasterize, kwargs...)
    for i ∈ xind
        x = zeros(length(ys)) .+ xs[i]
        lines!(ax, x, ys, zs[i, :]; color = :whitesmoke, linewidth)
    end 
    for i ∈ yind
        y = zeros(length(xs)) .+ ys[i]
        lines!(ax, xs, y, zs[:, i]; color = :whitesmoke, linewidth)
    end 
end

function surfacegrid_samples(vals, nsamples) 
    valslength = length(vals)
    inds = [ round(Int, (i - 1) * (valslength - 1) / (nsamples - 1) + 1) for i in 1:nsamples ]
    return inds
end 

function plotequilibriumsurface!(gl::GridLayout, betas, phis, Is; kwargs...)
    ax = Axis3(gl[1, 1]; azimuth = .7π)
    plotequilibriumsurface!(ax, betas, phis, Is; kwargs...)
    return (gl, ax)
end 

function plotequilibriumsurface!(ax::Axis3, betas, phis, Is; kwargs...)
    surfacegrid!(ax, betas, phis, Is; kwargs...)
end 

function labelequilibriumsurface!(ax::Axis3, labels::Vector{<:AbstractString}; 
        xlaboffset = 20, ylaboffset = 20, zlaboffset = 70, zlabrotation = 3π / 2
    )
    ax.xlabel = labels[1]
    ax.ylabel = labels[2]
    ax.zlabel = labels[3]
    ax.xlabeloffset = xlaboffset
    ax.ylabeloffset = ylaboffset
    ax.zlabeloffset = zlaboffset
    ax.zlabelrotation = zlabrotation
    formataxis!(ax)
end
