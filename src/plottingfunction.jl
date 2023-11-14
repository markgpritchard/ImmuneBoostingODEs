
# Makie is plotting y axes that are not completely vertical. A pull request is pending 
# to fix this. In the meantime, fix manually in this code 

function create_linepoints(
        pos_ext_hor,
        flipped::Bool, spine_width::Number, trimspine::Union{Bool, Tuple{Bool, Bool}}, tickpositions::Vector{Point2f}, tickwidth::Number)

    (position::Float32, extents::NTuple{2, Float32}, horizontal::Bool) = pos_ext_hor

    if trimspine isa Bool
        trimspine = (trimspine, trimspine)
    end

    if trimspine == (false, false) || length(tickpositions) < 2
        if horizontal
            y = position
            p1 = Point2f(extents[1] - 0.5spine_width, y)
            p2 = Point2f(extents[2] + 0.5spine_width, y)
            return [p1, p2]
        else
            x = position
            p1 = Point2f(x, extents[1] - 0.5spine_width)
            p2 = Point2f(x, extents[2] + 0.5spine_width)
            return [p1, p2]
        end
    else
        extents_oriented = last(tickpositions) > first(tickpositions) ? extents : reverse(extents)
        if horizontal
            y = position
            pstart = Point2f(-0.5f0 * tickwidth, 0)
            pend = Point2f(0.5f0 * tickwidth, 0)
            from = trimspine[1] ? tickpositions[1] .+ pstart : Point2f(extents_oriented[1] - 0.5spine_width, y)
            to = trimspine[2] ? tickpositions[end] .+ pend : Point2f(extents_oriented[2] + 0.5spine_width, y)
            return [from, to]
        else
            x = position
            pstart = Point2f(0, -0.5f0 * tickwidth)
            pend = Point2f(0, 0.5f0 * tickwidth)
            from = trimspine[1] ? tickpositions[1] .+ pstart : Point2f(x, extents_oriented[1] - 0.5spine_width)
            to = trimspine[2] ? tickpositions[end] .+ pend : Point2f(x, extents_oriented[2] + 0.5spine_width)
            return [from, to]
        end
    end

end
