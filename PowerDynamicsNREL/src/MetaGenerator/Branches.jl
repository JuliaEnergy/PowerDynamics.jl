"""
    get_pd_branch(arc<:branch, branches::Vector{Branch})

Add methods for each supported subtype of `Branch`. Returns an instance
of `AbstractLine`. There are potentially more than one branch objects
"""
function get_pd_branch end

function get_pd_branch(arc::Arc, lines::AbstractVector{Line})
    from = get_name(arc.from)
    to = get_name(arc.to)


    PiModelLine(; from, to, y, y_shunt_km, y_shunt_mk)
end
