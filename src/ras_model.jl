
function sirrrs!(du, u, p, t)
    # columns are age-groups, rows are compartments

    # if solely dividing by N_i, then when there is only one infant, if they happen to be 
    # infectious then they have an infeasible force of infection on the whole population. 
    # Therefore set a minimum N of one month's births. Will want to modify this slightly to 
    # make it continuously differentiable 
    λ = [ 
        p.τ * (
            p.c[i, 1] * sum(@view u[2, 1:5]) / (p.ν + sum(@view u[1:5, 1:5])) + 
            p.c[i, 2] * sum(@view u[2, 6:10]) / sum(@view u[1:5, 6:10]) + 
            p.c[i, 3] * sum(@view u[2, 11:15]) / sum(@view u[1:5, 11:15]) + 
            sum([ p.c[i, j] * u[2, (j + 12)] / sum(@view u[1:5, (j + 12)]) for j ∈ 4:16 ])
        ) * 365.25
        for i ∈ 1:16
    ]

    # age 0
    du[1, 1] = p.ν * sum(@view u[1:2, 17:21]) / sum(@view u[1:5, 17:21]) +              # S
        3 * p.ω * u[5, 1] - 
        λ[1] * u[1, 1]                                 
    du[2, 1] = λ[1] * u[1, 1] - p.γ * u[2, 1]                                           # I
    du[3, 1] = p.ν * sum(@view u[3, 17:21]) / sum(@view u[1:5, 17:21]) +                # R1
        p.γ * u[2, 1] + 
        λ[1] * p.ψ * (u[4, 1] + u[5, 1]) - 
        (3 * p.ω) * u[3, 1]   
    du[4, 1] = p.ν * sum(@view u[4, 17:21]) / sum(@view u[1:5, 17:21]) +                # R2
        3 * p.ω * u[3, 1] - 
        (3 * p.ω + λ[1] * p.ψ) * u[4, 1]                     
    du[5, 1] = p.ν * sum(@view u[5, 17:21]) / sum(@view u[1:5, 17:21]) +                # R3
        3 * p.ω * u[4, 1] - 
        (3 * p.ω + λ[1] * p.ψ) * u[5, 1]                     
    du[6, 1] = λ[1] * u[1, 1]                                                           # cumulative cases

    for i ∈ 2:15  # ages 1 to 14 
        if i <= 5 
            j = 1 
        elseif i <= 10 
            j = 2 
        else 
            j = 3 
        end
        du[1, i] = 3 * p.ω * u[5, i] - λ[j] * u[1, i]                                   # S
        du[2, i] = λ[j] * u[1, i] - p.γ * u[2, i]                                       # I
        du[3, i] = p.γ * u[2, i] +                                                      # R1
            λ[j] * p.ψ * (u[4, i] + u[5, i]) - 
            (3 * p.ω) * u[3, i]  
        du[4, i] = 3 * p.ω * u[3, i] - (3 * p.ω + λ[j] * p.ψ) * u[4, i]                 # R2
        du[5, i] = 3 * p.ω * u[4, i] - (3 * p.ω + λ[j] * p.ψ) * u[5, i]                 # R3
        du[6, i] = λ[j] * u[1, i]                                                       # cumulative cases
    end

    # ages 15 to 19 
    du[1, 16] = 3 * p.ω * u[5, 16] - λ[4] * u[1, 16] - (0.2 + p.μ[16]) * u[1, 16]       # S
    du[2, 16] = λ[4] * u[1, 16] - p.γ * u[2, 16] - (0.2 + p.μ[16]) * u[2, 16]           # I
    du[3, 16] = p.γ * u[2, 16] +                                                        # R1
        λ[4] * p.ψ * (u[4, 16] + u[5, 16]) - 
        (3 * p.ω) * u[3, 16] - 
        (0.2 + p.μ[16]) * u[3, 16]                                                      
    du[4, 16] = 3 * p.ω * u[3, 16] -                                                    # R2
        (3 * p.ω + λ[4] * p.ψ) * u[4, 16] - 
        (0.2 + p.μ[16]) * u[4, 16]             
    du[5, 16] = 3 * p.ω * u[4, 16] -                                                    # R3
        (3 * p.ω + λ[4] * p.ψ) * u[5, 16] - 
        (0.2 + p.μ[16]) * u[5, 16]             
    du[6, 16] = λ[4] * u[1, 16]                                                         # cumulative cases

    for i ∈ 17:27  # ages 20 to 74 
        j = i - 12
        du[1, i] = 3 * p.ω * u[5, i] -                                                  # S
            λ[j] * u[1, i] + 
            0.2 * u[1, i-1] - 
            (0.2 + p.μ[i]) * u[1, i]                 
        du[2, i] = λ[j] * u[1, i] -                                                     # I
            p.γ * u[2, i] + 
            0.2 * u[2, i-1] - 
            (0.2 + p.μ[i]) * u[2, i]                                  
        du[3, i] = p.γ * u[2, i] +                                                      # R1
            λ[j] * p.ψ * (u[4, i] + u[5, i]) - 
            (3 * p.ω) * u[3, i] + 
            0.2 * u[3, i-1] - 
            (0.2 + p.μ[i]) * u[3, i]      
        du[4, i] = 3 * p.ω * u[3, i] -                                                  # R2
            (3 * p.ω + λ[j] * p.ψ) * u[4, i] + 
            0.2 * u[4, i-1] - 
            (0.2 + p.μ[i]) * u[4, i]             
        du[5, i] = 3 * p.ω * u[4, i] -                                                  # R3
            (3 * p.ω + λ[j] * p.ψ) * u[5, i] + 
            0.2 * u[5, i-1] - 
            (0.2 + p.μ[i]) * u[5, i]             
        du[6, i] = λ[j] * u[1, i]                                                       # cumulative cases
    end

    # age >= 75 
    du[1, 28] = 3 * p.ω * u[5, 28] -                                                    # S
        λ[16] * u[1, 28] + 
        0.2 * u[1, 27] - 
        p.μ[28] * u[1, 28]                
    du[2, 28] = λ[16] * u[1, 28] -                                                      # I
        p.γ * u[2, 28] + 
        0.2 * u[2, 27] - 
        p.μ[28] * u[2, 28]                                  
    du[3, 28] = p.γ * u[2, 28] +                                                        # R1
        λ[16] * p.ψ * (u[4, 28] + u[5, 28]) - 
        (3 * p.ω) * u[3, 28] + 0.2 * u[3, 27] - 
        p.μ[28] * u[3, 28]      
    du[4, 28] = 3 * p.ω * u[3, 28] -                                                    # R2
        (3 * p.ω + λ[16] * p.ψ) * u[4, 28] + 
        0.2 * u[4, 27] - 
        p.μ[28] * u[4, 28]             
    du[5, 28] = 3 * p.ω * u[4, 28] -                                                    # R3
        (3 * p.ω + λ[16] * p.ψ) * u[5, 28] + 
        0.2 * u[5, 27] - 
        p.μ[28] * u[5, 28]             
    du[6, 28] = λ[16] * u[1, 28]                                                        # cumulative cases
end

function makecontactmatrix(source)
    contacts_df = CSV.read(source, DataFrame; header=false)
    return Matrix(contacts_df)
end

function makeoutputmatrices(sol)
    t = sol.t 
    susceptible = makeoutputmatrix_susceptible(sol)
    infectious = makeoutputmatrix_infectious(sol)
    resistant = makeoutputmatrix_resistant(sol)
    N = makeoutputmatrix_N(sol)
    cumcases = makeoutputmatrix_cumcases(sol)
    incidentcases = makeoutputmatrix_incidentcases(sol)
    return @ntuple t susceptible infectious resistant N cumcases incidentcases
end

function _makeoutputmatrix(sol::SciMLBase.AbstractODESolution{T, <:Any, <:Any}, ind) where T
    mat = Matrix{T}(undef, length(sol), 7)
    _fillputmatrix!(mat, sol, ind)
    return mat
end

function _fillputmatrix!(mat, sol, ind::Integer)
    for i ∈ eachindex(sol) 
        mat[i, 1] = sol[i][ind, 1]
        mat[i, 2] = sum(@view sol[i][ind, 2:5])
        mat[i, 3] = sum(@view sol[i][ind, 6:15])
        mat[i, 4] = sum(@view sol[i][ind, 16:21])
        mat[i, 5] = sum(@view sol[i][ind, 22:25])
        mat[i, 6] = sum(@view sol[i][ind, 26:27])
        mat[i, 7] = sol[i][ind, 28]
    end
end

function _fillputmatrix!(mat, sol, inds)
    for i ∈ eachindex(sol) 
        mat[i, 1] = sum(@view sol[i][inds, 1])
        mat[i, 2] = sum(@view sol[i][inds, 2:5])
        mat[i, 3] = sum(@view sol[i][inds, 6:15])
        mat[i, 4] = sum(@view sol[i][inds, 16:21])
        mat[i, 5] = sum(@view sol[i][inds, 22:25])
        mat[i, 6] = sum(@view sol[i][inds, 26:27])
        mat[i, 7] = sum(@view sol[i][inds, 28])
    end
end

makeoutputmatrix_susceptible(sol) = _makeoutputmatrix(sol, 1)
makeoutputmatrix_infectious(sol) = _makeoutputmatrix(sol, 2)
makeoutputmatrix_resistant(sol) = _makeoutputmatrix(sol, 3:5)
makeoutputmatrix_N(sol) = _makeoutputmatrix(sol, 1:5)
makeoutputmatrix_cumcases(sol) = _makeoutputmatrix(sol, 6)

function makeoutputmatrix_incidentcases(sol)
    incidentcases = makeoutputmatrix_cumcases(sol)
    ℓ = size(incidentcases, 1)
    for i ∈ axes(incidentcases, 1), j ∈ axes(incidentcases, 2)
        k = ℓ - i + 1
        if k == 1 
            incidentcases[k, j] = zero(incidentcases[k, j])
        else
            incidentcases[k, j] += -incidentcases[(k - 1), j]
        end
    end
    return incidentcases
end

memosolver(prob, alg; kwargs...) = @memoize solve(prob, alg; kwargs...)

function solvesamples(prob, samplevalues; kwargs...)
    nsamples = size(samplevalues, 1)
    sample1 = solvesample(prob, samplevalues, 1; kwargs...)
    output = Array{Float64}(undef, size(sample1)..., nsamples)
    output[:, :, 1] .= sample1
    for i ∈ 2:nsamples
        output[:, :, i] .= solvesample(prob, samplevalues, i; kwargs...)
    end
    means = Array{Float64}(undef, size(sample1)...)
    lq = Array{Float64}(undef, size(sample1)...)
    uq = Array{Float64}(undef, size(sample1)...)
    for i ∈ axes(output, 1), j ∈ axes(output, 2)
        means[i, j] = mean(output[i, j, :])
        lqv, uqv = quantile(output[i, j, :], [ 0.05, 0.95 ])
        lq[i, j] = lqv
        uq[i, j] = uqv
    end
    return @ntuple means lq uq
end

function solvesample(prob, samplevalues, ind::Integer; kwargs...)
    τ = samplevalues.τ[ind]
    ψ = samplevalues.ψ[ind]
    ϕ = samplevalues.ϕ[ind]
    return solvesample(prob, τ, ψ, ϕ; kwargs...)
end

function solvesample(prob, τ, ψ, ϕ;
    births, callbackset, contacts_alllocations, mortality, saveat, u0, 
    abstol=1e-15, alg=Vern9(lazy=false), gamma=48.7, maxiters=5e4, omega=0.913
)
    p = SirrrsParameters(τ, contacts_alllocations, gamma, births[1], mortality, ψ, omega)
    sol = memosolver(prob, alg; p, u0, callback=callbackset, saveat, abstol, maxiters,)
    incidentcases = makeoutputmatrix_incidentcases(sol)
    return incidentcases .* ϕ
end
