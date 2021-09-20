"""
    optimize_OptStoic(os_model, return_solution = true)

Optimizes an previously built optstoic model. Returns optimized model `os_model` and a
dictionary with the optimized stoichiometry `res_dc` if `return_solution` is `true`.
"""
# function optimize_OptStoic(os_model, return_solution = true)
# optimize!(os_model)
# if return_solution == false
#     return os_model
# else
#     res_dc = Dict()
#     res_dc[name(os_model[:S])] = value(os_model[:S])
#     for i in os_model[:T]
#         res_dc[name(i)[1:end-3]] = value(i)
#     end
#     for i in os_model[:CoR]
#         res_dc[name(i)[1:end-3]] = value(i)
#     end
# end
# return os_model, res_dc
# end

"""
    MinFlux procedure to optimize model and returns optimized model. Currently work in progress, model build needs to be fixed.
"""
function MinFlux(mf_model)
    optimize!(mf_model)
    return mf_model
end
