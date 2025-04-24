
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Effect of boosting on equilibria 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Compartments over time 

function plotequilibriumc!(fig::Figure, basic5, basic1_5, basic5_; row=2)
    gc = GridLayout(fig[row, 1:3])
    plotequilibriumc!(gc, basic5, basic1_5, basic5_)
    labelplots!("C", gc; rows = 1) 
end 

function plotequilibriumc!(gc::GridLayout, basic5, basic1_5, basic5_)
    yt = [ 0.4:0.2:0.8, [ 0, 0.02 ], 0:0.2:0.4 ]

    Label(gc[0, 1], L"$\psi=0.5$"; halign=:left, fontsize=11.84, tellwidth=false)
    axs1 = [ Axis(gc[i, 1]; yticks=yt[i]) for i ∈ 1:3 ]
    lines!(axs1[1], basic5[:gt], basic5[:S]; color=COLOUR_S, linewidth=1,)
    lines!(axs1[2], basic5[:gt], basic5[:I]; color=COLOUR_I, linewidth=1,)
    lines!(axs1[3], basic5[:gt], basic5[:R1]; color=COLOURVECTOR[3], linewidth=1,)
    lines!(axs1[3], basic5[:gt], basic5[:R2]; color=COLOURVECTOR[4], linewidth=1,)
    lines!(axs1[3], basic5[:gt], basic5[:R3]; color=COLOURVECTOR[5], linewidth=1,)
    
    Label(gc[0, 2], L"$\psi=1.5$"; halign=:left, fontsize=11.84, tellwidth=false)
    axs2 = [ Axis(gc[i, 2]; yticks=yt[i]) for i ∈ 1:3 ]
    lines!(axs2[1], basic1_5[:gt], basic1_5[:S]; color=COLOUR_S, linewidth=1,)
    lines!(axs2[2], basic1_5[:gt], basic1_5[:I]; color=COLOUR_I, linewidth=1,)
    lines!(axs2[3], basic1_5[:gt], basic1_5[:R1]; color=COLOURVECTOR[3], linewidth=1,)
    lines!(axs2[3], basic1_5[:gt], basic1_5[:R2]; color=COLOURVECTOR[4], linewidth=1,)
    lines!(axs2[3], basic1_5[:gt], basic1_5[:R3]; color=COLOURVECTOR[5], linewidth=1,)
    
    Label(gc[0, 3], L"$\psi=5$"; halign=:left, fontsize=11.84, tellwidth=false)
    axs3 = [ Axis(gc[i, 3]; yticks=yt[i]) for i ∈ 1:3 ]
    lines!(axs3[1], basic5_[:gt], basic5_[:S]; color=COLOUR_S, linewidth=1,)
    lines!(axs3[2], basic5_[:gt], basic5_[:I]; color=COLOUR_I, linewidth=1,)
    lines!(axs3[3], basic5_[:gt], basic5_[:R1]; color=COLOURVECTOR[3], linewidth=1,)
    lines!(axs3[3], basic5_[:gt], basic5_[:R2]; color=COLOURVECTOR[4], linewidth=1,)
    lines!(axs3[3], basic5_[:gt], basic5_[:R3]; color=COLOURVECTOR[5], linewidth=1,)
    
    leg = Legend(
        gc[3, 4], 
        [ LineElement(; color=c, linewidth=1,) for c ∈ COLOURVECTOR[3:5] ], 
        [ "$ℓ" for ℓ ∈ [ L"$R_1$", L"$R_2$", L"$R_3$" ]]; 
        titlefont=:regular, orientation=:vertical,
    )

    linkxaxes!(axs1..., axs2..., axs3...)
    for j ∈ 1:3 
        linkyaxes!(axs1[j], axs2[j], axs3[j])
        for (i, axs) ∈ enumerate([ axs1, axs2, axs3 ])
            formataxis!(
                axs[j]; 
                hidex=(j != 3), hidexticks=(j != 3), hidey=(i != 1), hideyticks=(i != 1), 
                trimspines=true,
            )
            hidespines!(axs[j], :t)
            hidespines!(axs[j], :r)
            if i != 1 hidespines!(axs[j], :l) end
            if j == 1 setvalue!(axs[j], 0, 0.4) end
            if j != 3 hidespines!(axs[j], :b) end
        end 
    end 
    formataxis!(leg; horizontal = false)

    Label(gc[4, 1:3], "Time, years"; fontsize = 11.84, tellwidth = false)
    Label(gc[1, 0], "Susceptible"; fontsize = 11.84, tellheight = false)
    Label(gc[2, 0], "Infectious"; fontsize = 11.84, tellheight = false)
    Label(gc[3, 0], "Resistant"; fontsize = 11.84, tellheight = false)
    
    for c ∈ [ 1, 4 ] colgap!(gc, c, 5) end
    for r ∈ [ 1, 4 ] rowgap!(gc, r, 5) end
    for r ∈ [ 2, 3 ] rowgap!(gc, r, 7) end
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fourier transforms 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotfourier!(
    fig, 
    unforced6mfreqs, 
    unforced365dfreqs,
    unforced400dfreqs, 
    unforced25freqs, 
    unforced6mdensities, 
    unforced365ddensities,
    unforced400ddensities, 
    unforced25densities, 
    forced6mfreqs, 
    forced365dfreqs,
    forced400dfreqs, 
    forced25freqs, 
    forced6mdensities, 
    forced365ddensities,
    forced400ddensities, 
    forced25densities, 
    simparms
)
    gl = GridLayout(fig[1, 1])

    Label(gl[-1, 1], "Immune duration:"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(gl[0, 1], "6 months"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(gl[0, 2], "1.1 years"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(gl[0, 3], "2.5 years"; fontsize=11.84, halign=:left, tellwidth=false)

    _plotfouriera!(
        gl, 
        unforced6mfreqs, 
        unforced365dfreqs,
        unforced400dfreqs, 
        unforced25freqs, 
        unforced6mdensities, 
        unforced365ddensities,
        unforced400ddensities, 
        unforced25densities, 
        simparms
    )
    _plotfourierb!(
        gl, 
        forced6mfreqs, 
        forced365dfreqs,
        forced400dfreqs, 
        forced25freqs, 
        forced6mdensities, 
        forced365ddensities,
        forced400ddensities, 
        forced25densities, 
        simparms
    ) 

    labelplots!([ "A", "B"], gl; rows=[ -1, 2 ])
    for r ∈ [ 1, 2, 4 ] rowgap!(gl, r, 5) end
    for c ∈ [ 1, 5 ] colgap!(gl, c, 5) end
end

function _plotfouriera!(
    gl, 
    unforced6mfreqs, 
    unforced365dfreqs,
    unforced400dfreqs, 
    unforced25freqs, 
    unforced6mdensities, 
    unforced365ddensities,
    unforced400ddensities, 
    unforced25densities, 
    simparms
)
    @unpack psis = simparms
    freqvector = [ unforced6mfreqs, unforced365dfreqs, unforced400dfreqs, unforced25freqs ]
    densityvector = [ unforced6mdensities, unforced365ddensities, unforced400ddensities, unforced25densities ]

    axs = [ Axis(gl[1, i]; yscale=log) for i ∈ 1:4 ]
    i = 1
    for (ax, freqs, densities) ∈ zip(axs, freqvector, densityvector)
        inds = findall(x -> .1 <= x <= 365 / 21, freqs) 
        hm = heatmap!(
            ax, psis, freqs[inds], log10.(densities[inds, :]'); 
            colormap=:CMRmap, colorrange=( -6, 0 ), rasterize=10,
        )
        ax.yticks = ( 
            [ 0.2, 0.5, 1, 2, 4, 365 / 35 ], 
            [ "5 y", "2 y", "1 y", "6 m", "3 m", "5 w" ] 
        )
        formataxis!(ax; hidex=true, hidey=(i != 1))
        if i == 1 
            co = Colorbar(gl[1:2, 4], hm) 
            formataxis!(co)
        end
        i += 1
    end
    for i ∈ 1:2
        Label(gl[i, 0], "Period"; fontsize=11.84, rotation=π/2, tellheight=false)
    end
    Label(
        gl[1:2, 5], L"$\log_{10}$ spectral density"; 
        fontsize=11.84, rotation=-π/2, tellheight=false
    )
end

function _plotfourierb!(
    gl, 
    forced6mfreqs, 
    forced365dfreqs,
    forced400dfreqs, 
    forced25freqs, 
    forced6mdensities, 
    forced365ddensities,
    forced400ddensities, 
    forced25densities, 
    simparms
) 
    @unpack psis = simparms
    freqvector = [ forced6mfreqs, forced365dfreqs, forced400dfreqs, forced25freqs ]
    densityvector = [ forced6mdensities, forced365ddensities, forced400ddensities, forced25densities ]

    axs = [ Axis(gl[2, i]; yscale=log) for i ∈ 1:4 ]
    i = 1
    for (ax, freqs, densities) ∈ zip(axs, freqvector, densityvector)
        inds = findall(x -> 0.1 <= x <= 365 / 21, freqs) 
        hm = heatmap!(
            ax, psis, freqs[inds], log10.(densities[inds, :]'); 
            colormap = :CMRmap, colorrange = ( -6, 0 )
        )
        ax.yticks = (
            [ 0.2, 0.5, 1, 2, 4, 365 / 35 ], 
            [ "5 y", "2 y", "1 y", "6 m", "3 m", "5 w" ]
        )
        formataxis!(ax; hidey=(i != 1))
        i += 1
    end
    Label(
        gl[3, 1:3], L"Boosting coefficient, $\psi$"; 
        fontsize=11.84, tellwidth=false
    )
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Non-pharmaceutical interventions plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotnpi!(fig, npisim_phi0, npisim_phi5, npisim_phi13_2, npiparms) 
    @unpack reductiontime = npiparms

    ga = GridLayout(fig[1, 1])
    axs = [ Axis(ga[i, j]) for i ∈ [ 1, 3, 5 ], j ∈ [ 1, 3 ] ]
    jax1 = Axis(ga[1:5, 1])
    jax2 = Axis(ga[1:5, 3])
    for (i, m) ∈ enumerate([ npisim_phi0, npisim_phi5, npisim_phi13_2 ])
        @unpack cases, compartments = m 
        lines!(
            axs[i, 1], compartments[:gt][2:end], cases; 
            color=COLOUR_I, linewidth=1,
        )
        lines!(
            axs[i, 2], compartments[:gt], compartments[:S]; 
            color=COLOUR_S, linewidth=1, label=L"$S$",
        )
        lines!(
            axs[i, 2], compartments[:gt], compartments[:R1]; 
            color=COLOURVECTOR[3], linewidth=1, label=L"$R_1$",
        )
        lines!(
            axs[i, 2], compartments[:gt], compartments[:R2]; 
            color=COLOURVECTOR[5], linewidth=1, label=L"$R_2$",
        )
        lines!(
            axs[i, 2], compartments[:gt], compartments[:R3]; 
            color=COLOURVECTOR[6], linewidth=1, label=L"$R_3$",
        )

    end
    for ax ∈ [ jax1, jax2 ]
        for x ∈ 0:1:10  
            vlines!(
                ax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
            )
        end
        vspan!(ax, reductiontime, reductiontime + 1, color=( :gray, 0.1 )) 
        formataxis!(
            ax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b ),
        )
    end
    leg = Legend(ga[7, 3], axs[2, 2]; padding = ( 5, 5, 2, 2 ))

    linkxaxes!(axs..., jax1, jax2)
    for j ∈ 1:2 linkyaxes!(axs[:, j]...) end
    for i ∈ 1:3, j ∈ 1:2
        formataxis!(axs[i, j]; hidex=(i != 3), hidexticks=(i != 3), trimspines=true)
        hidespines!(axs[i, j], :t)
        hidespines!(axs[i, j], :r)
        i != 3 && hidespines!(axs[i, j], :b)
    end
    formataxis!(leg)

    Label(
        ga[0, 1], L"Mean $\mathcal{R}_0=1.215$, with seasonal forcing $± 10%,\ \psi=0$"; 
        fontsize=11.84, halign=:left, tellwidth=false
    )
    Label(
        ga[2, 1], L"Mean $\mathcal{R}_0=1.285$, with seasonal forcing $± 8.2%,\ \psi= 5$";
        fontsize=11.84, halign=:left, tellwidth=false
    )
    Label(
        ga[4, 1], L"$\mathcal{R}_0=1.6$, with no seasonal forcing, $\psi= 13.2$";
        fontsize=11.84, halign=:left, tellwidth=false
    )

    for i ∈ [ 1, 3 ]
        Label(ga[6, i], "Time, years"; fontsize=11.84, tellwidth=false)
    end
    for i ∈ [ 1, 3, 5 ]
        Label(ga[i, 0], "Weekly incidence"; fontsize=11.84, rotation=π/2, tellheight=false)
        Label(ga[i, 2], "Proportion"; fontsize=11.84, rotation=π/2, tellheight=false)
    end
    for r ∈ [ 1, 3, 5, 6, 7 ] rowgap!(ga, r, 5) end
    for c ∈ [ 1, 3 ] colgap!(ga, c, 5) end
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-stratified respiratory syncytial virus data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotrsvage!(fig, agedata)
    ms = [ 0, 91, 183, 275, 365 ] ./ 365
    colours = [ ( :gray, 0.4 ); [ COLOURVECTOR[i] for i ∈ 1:3 ] ]
    agelabels = Dict(
        "<1 years" => L"$<1$ years",
        "1-4 years" => L"$1$–$4$ years",
        "5-14 years" => L"$5$–$14$ years",
        "15-44 years" => L"$15$–$44$ years",
        "45-64 years" => L"$45$–$64$ years",
        "65-74 years" => L"$65$–$74$ years",
        "75+ years" => L"$\ge 75$ years",
    )
    xticks = ( ms, [ "April", "July", "Oct", "Jan", "April" ] )

    axs = [ i <= 4 ? Axis(fig[i, 1]; xticks) : Axis(fig[i-4, 2]; xticks) for i ∈ 1:7 ]
    axv1 = Axis(fig[1:4, 1]; xticks)
    axv2 = Axis(fig[1:3, 2]; xticks)
    for (i, y) ∈ enumerate(collect(2016:1:2022))
        d = subset(agedata, :Year => x -> x.== y)
        for (j, a) ∈ enumerate(unique(agedata.AgeGroup))
            d2 = subset(d, :AgeGroup => x -> x.== a )
            if y <= 2019 
                color = colours[1] 
            else 
                color = colours[i-3] 
            end
            lines!(
                axs[j], d2.FractionDate, getproperty(d2, Symbol("Rate$a")); 
                color, label="$y", linewidth=1,
            ) 
            if y == 2022 
                maxy = maximum(getproperty(d2, Symbol("Rate$a")))
                text!(
                    axs[j], 0, 0.9 * maxy; 
                    text="$(agelabels[a])", fontsize=11.84, align=( :left, :top )
                )
            end
        end
    end       
    
    for ax ∈ [ axv1, axv2 ], x ∈ ms
        vlines!(
            ax, x; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
        )
    end

    for i ∈ 1:7
        formataxis!(
            axs[i]; 
            hidex=(i != 4 && i != 7), hidexticks=(i != 4 && i != 7), trimspines=true
        ) 
        hidespines!(axs[i], :t)
        hidespines!(axs[i], :r)
        if i != 4 && i != 7 
            hidespines!(axs[i], :b) 
        end
    end

    for ax ∈ [ axv1, axv2 ]
        formataxis!(
            ax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        ) 
    end

    linkxaxes!(axs..., axv1, axv2)

    # create legend manually 
    lineelements = [ LineElement(; color=c) for c ∈ colours ]
    labels = [ "2016-19", "2020", "2021", "2022" ]
    leg = Legend(
        fig[4:5, 2], lineelements, labels, "Year:"; 
        titlefont=:regular, tellheight=false, tellwidth=false, nbanks=2, valign=:bottom,
    )
    leg.framevisible = false
    leg.labelsize = 10
    leg.titlesize = 12
    leg.titleposition = :top

    Label(
        fig[5, 1], "Month"; 
        fontsize=11.84, tellwidth=false,
    )
    Label(
        fig[4, 2], "Month"; 
        fontsize=11.84, tellheight=false, tellwidth=false, valign=:top,
    )
    Label(
        fig[1:4, 0], "Annual cumulative incidence"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    _formatplotrsvage!(fig)
end

_formatplotrsvage!(fig::Figure) = __formatplotrsvage!(fig.layout)
_formatplotrsvage!(gl::GridLayout) = __formatplotrsvage!(gl)

function __formatplotrsvage!(fl)
    colgap!(fl, 1, 5) 
    rowgap!(fl, 3, 5) 
    rowgap!(fl, 4, 5) 
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Susceptible--infectious planes 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotsi!(
    fig::Figure, 
    stablesim, 
    unstablein12, 
    unstableout12,
    unstablein5, 
    unstableout5, 
    limitcycle12, 
    limitcycle5;
    row=1
)
    _plotsi!(
        fig, 
        stablesim, 
        unstablein12, 
        unstableout12, 
        unstablein5, 
        unstableout5, 
        limitcycle12, 
        limitcycle5;
        row
    )
    colgap!(fig.layout, 1, 5) 
    for r ∈ row:(row + 1) rowgap!(fig.layout, r, 5) end
end

function plotsi!(
    fig::GridLayout, 
    stablesim, 
    unstablein12,
    unstableout12, 
    unstablein5, 
    unstableout5, 
    limitcycle12,
    limitcycle5; 
    row=1
)
    _plotsi!(
        fig, 
        stablesim, 
        unstablein12,
        unstableout12, 
        unstablein5, 
        unstableout5, 
        limitcycle12, 
        limitcycle5; 
        row
    )
    colgap!(fig, 1, 5) 
    rowgap!(fig, row + 1, 5) 
end

function _plotsi!(
    fig, 
    stablesim, 
    unstablein12, 
    unstableout12, 
    unstablein5, 
    unstableout5, 
    limitcycle12, 
    limitcycle5;
    row=1
)
    ax1 = Axis(fig[row, 1])
    lines!(ax1, stablesim[:S], stablesim[:I]; color=COLOURVECTOR[1], linewidth=1,)
    for ℓ ∈ [ 3050, 5900 ]
        x = stablesim[:S][ℓ]
        y = stablesim[:I][ℓ]
        u = x - stablesim[:S][ℓ-1] 
        v = y - stablesim[:I][ℓ-1] 
        arrows!(ax1, [x], [y], [u], [v]; color=COLOURVECTOR[1])
    end

    ax2 = Axis(fig[row, 2])
    lines!(ax2, limitcycle12[:S], limitcycle12[:I]; color=:black, linewidth=1,)
    for (i, r) ∈ enumerate([ unstablein12, unstableout12 ])
        lines!(ax2, r[:S], r[:I]; color=COLOURVECTOR[i+1])
        ℓ = length(r[:S])
        x = r[:S][ℓ]
        y = r[:I][ℓ]
        u = x - r[:S][ℓ-1] 
        v = y - r[:I][ℓ-1] 
        arrows!(ax2, [x], [y], [u], [v]; color=COLOURVECTOR[i+1])
    end

    ax3 = Axis(fig[row, 3])
    lines!(ax3, limitcycle5[:S], limitcycle5[:I]; color=:black, linewidth=1,)
    for (i, r) ∈ enumerate([ unstablein5, unstableout5 ])
        lines!(ax3, r[:S], r[:I]; color=COLOURVECTOR[i+1])
        ℓ = length(r[:S])
        x = r[:S][ℓ]
        y = r[:I][ℓ]
        u = x - r[:S][ℓ-1] 
        v = y - r[:I][ℓ-1] 
        arrows!(ax3, [x], [y], [u], [v]; color=COLOURVECTOR[i+1])
    end

    linkaxes!(ax1, ax2, ax3)

    formataxis!(
        ax1; 
        hidespines=( :r, :t ), trimspines=true,
    )
    formataxis!(
        ax2; 
        hidey=true, hideyticks=true, hidespines=(:l, :r, :t), trimspines=true,
    )
    formataxis!(
        ax3; 
        hidey=true, hideyticks=true, hidespines=(:l, :r, :t), trimspines=true,
    )
    Label(fig[(row + 1), 1:3], "Susceptible"; fontsize=11.84, tellwidth=false)
    Label(fig[row, 0], "Infectious"; fontsize=11.84, tellheight=false)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Strigency Index
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotstringency!(fig, crgtdata, reduceday, increaseday; kwargs...) 
    ax = Axis(fig[1, 1])
    vspan!(ax, reduceday, increaseday; color=( :gray, 0.1 ))
    lines!(ax, crgtdata.Date, crgtdata.StringencyIndex_Average; color=COLOURVECTOR[1])
    Label(fig[2, 1], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig[1, 0], "Stringency index"; fontsize=11.84, rotation=π/2, tellheight=false)
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    formataxis!(ax; kwargs...)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot MCMC trace plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotchains(data::DataFrame; size=( 400, 800 ), kwargs...)
    @unpack colnames, plotnames_ind = _processplotchains(data; kwargs...)
        
    fig = Figure(; size)
    ax = [ Axis(fig[i, 1]) for i ∈ eachindex(plotnames_ind) ]
    for (j, chainid) ∈ enumerate(unique(data.chain))
        inds = findall(x -> x == chainid, data.chain)
        for (i, k) ∈ enumerate(plotnames_ind) 
            lines!(ax[i], getproperty(data, colnames[k])[inds]; color=COLOURVECTOR[j])
            Label(fig.layout[i, 0], "$(colnames[k])"; rotation=π/2, tellheight=false)
        end
    end
    
    return fig
end

function _processplotchains(data; logdensity="log_density")
    # "log_density" is the label given by `Pigeons` output. Turing labels it "lp".
    colnames = names(data)
    lp_ind = findall(x -> x == logdensity, colnames)
    _plotnames_ind = findall(x -> x ∉ _NOPLOTNAMES, colnames)
    plotnames_ind = [ lp_ind; _plotnames_ind ]
    return @ntuple colnames plotnames_ind
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot simulations using MCMC values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plotfittedsimulationquantiles!(ax, data, plotvals, saveat)
    scatter!(
        ax, saveat[2:end], data.Cases; 
        color=:black, markersize=2,
    )
    lines!(
        ax, saveat[2:end], [ y[2] for y ∈ plotvals ]; 
        color=COLOUR_I, linewidth=1,
    )
    band!(
        ax, saveat[2:end], [ y[1] for y ∈ plotvals ], [ y[3] for y ∈ plotvals ]; 
        color=( COLOUR_I, 0.5 ),
    )
end
