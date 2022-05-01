function _narrow_round(v, thres)
    if v <= thres
        return 0
    else
        return v
    end
end

"""
    optimize_OptStoic(os_model, return_solution = true)

Optimises an previously built optstoic model. Returns optimised model `os_model` and a
dictionary with the optimised stoichiometry `res_dc`.
"""
function optimize_OptStoic(os_model)
    optimize!(os_model)
    return os_model
end

function solution_OptStoic(os_model)
    if is_solved(os_model)
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
    else
        return false
    end
end

"""
    DEV NOTE: Necessary for MinFlux/MinRxn procedure, currently work in progress. Function might change severely.
    MinFlux procedure to optimize model and returns optimized model. Currently work in progress, model build needs to be fixed.
"""
function MinFlux(mf_model, database, res = false)
    optimize!(mf_model)
    if res
        res_dc = Dict()
        for rxn in 1:length(reactions(database))
            if value(mf_model[:V][rxn]) != 0
                res_dc[reactions(database)[rxn]] = _narrow_round(value(mf_model[:V][rxn]), 50^-5)
            end
        end
        return mf_model, res_dc
    else
        return mf_model
    end
end
