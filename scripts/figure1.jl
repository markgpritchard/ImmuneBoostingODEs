
# Create the model flowchart (Figure 1 in the manuscript)

using DrWatson
@quickactivate :ImmuneBoostingODEs
using CairoMakie
CairoMakie.activate!() # allows figures to be saved as vector files

# Make the plot
introplot = Figure(; figure_padding = 1, resolution = ( 400, 200 ))
let 
    fontsize = 11.84 # font size for parameter labels on arrows 

    ga = GridLayout(introplot[1, 1])
    ax = Axis(ga[1, 1]) # axis that figure will be plotted on

    for i ∈ 1:5
        color = ( [ COLOUR_S, COLOUR_I, COLOUR_R, COLOUR_R, COLOUR_R ][i], .2 )
        label = [ "S", "I", "R₁", "R₂", "R₃" ][i]
        # generate coloured boxes for each compartment
        poly!(ax, Rect(i - .25, -.25, .5, .5); color)
        # label each box
        text!(ax, i, 0; text = label, align = ( :center, :center ))
    end 
    # arrows between boxes
    for x ∈ 2:4 arrows!(ax, [ x + .28 ], [ 0 ], [ .4 ], [ 0 ]; color = :black) end
    arrows!(ax, [ .52 ], [ .52 ], [ .2 ], [ -.2 ]; color = :black)
    for x ∈ 1:5 arrows!(ax, [ x + .28 ], [ -.28 ], [ .2 ], [ -.2 ]; color = :black) end
    # dashed lines
    xs1 = collect(1.28:.05:1.63); ys1 = zeros(length(xs1))
    linesegments!(ax, xs1, ys1; color = :red)
    arrows!(ax, [ 1.63 ], [ 0 ], [ .05 ], [ 0 ]; color = :red)
    xs2 = collect(4:-.05:3.05); ys2 = @. .52 - (xs2 - 3.5)^2
    linesegments!(ax, xs2, ys2; color = :red)
    xs3 = collect(5:-.05:3.05); ys3 = @. .57 - .3 * (xs3 - 4)^2
    linesegments!(ax, xs3, ys3; color = :red)
    arrows!(ax, [ 3.05 ], [ .2975 ], [ -.05 ], [ -.0275 ]; color = :red)
    lines!(ax, [ 5, 5 ], [ -.28, -.93 ]; color = :black)
    lines!(ax, [ 1, 5 ], [ -.93, -.93 ]; color = :black)
    arrows!(ax, [ 1 ], [ -.93 ], [ 0 ], [ .58 ]; color = :black)
    # label the arrows / lines
    text!(ax, 1.5, .1; text = "β(t) I(t)", color = :red, fontsize, align = ( :center, :bottom ))
    text!(ax, 2.5, .1; text = "γ", color = :black, fontsize, align = ( :center, :bottom ))
    for x ∈ 3:4
        text!(ax, x + .5, .1; 
            text = "3ω", color = :black, fontsize, align = ( :center, :bottom ))
    end
    text!(ax, .72, .45; text = "μ", color = :black, fontsize, align = ( :left, :bottom ))
    for x ∈ 1:5
        text!(ax, x + .48, -.35; text = "μ", color = :black, fontsize, align = ( :left, :bottom ))
    end
    text!(ax, 3.8, .65; text = "ψ β(t) I(t)", color = :red, fontsize, align = ( :center, :bottom ))
    text!(ax, 3, -1.03; text = "3ω", color = :black, fontsize, align = ( :center, :top ))

    formataxis!(ax; hidex = true, hidexticks = true, hidey = true, hideyticks = true, 
        hidespines = ( :l, :t, :r, :b ))
    setvalue!(ax, (.5, -1.1 )); setvalue!(ax, (5.5 , .75 ))
end 

# view the plot
introplot
# save the plot
safesave(plotsdir("introplot.pdf"), introplot)
