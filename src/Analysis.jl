"""
    optimize_OptStoic(os_model, return_solution = true)

Optimises an previously built optstoic model. Returns optimised model `os_model` and a
dictionary with the optimised stoichiometry `res_dc`.
"""
function optimize_OptStoic(os_model)
    optimize!(os_model)
    res_dc = Dict()
    res_dc[name(os_model[:S])] = value(os_model[:S])
    for i in os_model[:T]
        res_dc[name(i)[1:end-3]] = value(i)
    end
    for i in os_model[:CoR]
        res_dc[name(i)[1:end-3]] = value(i)
    end

    return os_model, res_dc
end

"""
    DEV NOTE: Necessary for MinFlux/MinRxn procedure, currently work in progress. Function might change severely.
    MinFlux procedure to optimize model and returns optimized model. Currently work in progress, model build needs to be fixed.
"""
function MinFlux(mf_model)
    optimize!(mf_model)
    return mf_model
end
