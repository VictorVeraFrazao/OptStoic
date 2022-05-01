"""
    optimize_OptStoic(os_model, return_solution = true)

Returns optimized OptStoic model.
"""
function optimize_OptStoic(os_model)
    optimize!(os_model)
    return os_model
end

"""
    solution_OptStoic(os_model)
Return dictionary with overall stoichiometry (metabolite => stoichiometric coefficient), if `os_model` is solved. If not it returns nothing.
"""
function solution_OptStoic(os_model)
    is_solved(os_model) || return nothing
    res_dc = Dict()
    for i in 1:length(os_model[:S])
        tstring = replace(name(os_model[:S][i]), string("[", i, "]") => "")
        res_dc[tstring] = value(os_model[:S][i])
    end
    for i in 1:length(os_model[:T])
        tstring = replace(name(os_model[:T][i]), string("[", i, "]") => "")
        res_dc[tstring] = value(os_model[:T][i])
    end
    for i in 1:length(os_model[:CoR])
        tstring = replace(name(os_model[:CoR][i]), string("[", i, "]") => "")
        res_dc[tstring] = value(os_model[:CoR][i])
    end
    return res_dc
end
"""
    MinFlux(mf_model)

Returns optimized MinFlux model.
"""
function MinFlux(mf_model)
    optimize!(mf_model)
    return mf_model
end

"""
    MinFlux_solution(mf_model)
Return dictionary with assembled pathway (reactions => flux), if `mf_model` is solved. If not it returns nothing. The basing metabolic model is passed via database.
"""
function MinFlux_solution(database, mf_model)
    is_solved(mf_model) || return nothing
    sol = Dict()
    for rxn in 1:length(reactions(database))
        if value(mf_model[:V][rxn]) != 0
            sol[reactions(database)[rxn]] = value(mf_model[:V][rxn])
        end
    end
    return sol
end