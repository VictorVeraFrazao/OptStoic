###DEPRECATED###VALID FUNCTIONS MOVE TO PREPROCESSING###WILL BE REMOVED SOON###




"""
    function adjust_model(database, variation_tres =.5, dGr_dict =Dict())

DEV NOTE: Necessary for MinFlux/MinRxn procedure, currently work in progress. Function might change severely.
Reduces model to usable reactions. It excludes exchange reactions, those where no Gibbs free energy could be calculated or with too large standard deviation, depending on the given coefficient of variation threshold (`variation_tres`). 
"""
function adjust_model(database, variation_tres = 0.9)
    dGr_dict, dGr_bounds = collect_dGr_bounds(database)

    for rxn in reactions(database)
        if !reaction_mass_balanced(database, rxn) #excluding exchanges
            database = remove_reaction(database, rxn)
        end
    end
    mod_bounds = OrderedDict()
    for (k, v) in dGr_bounds
        if k in reactions(database)
            mod_bounds[k] = v
        end
    end
    for rxn in keys(mod_bounds)
        if typeof(dGr_dict[rxn]) == String #excluding reactions without calculated Gibbs energies
            mod_bounds[rxn] = (0.0, 0.0)
        else
            try
                if dGr_dict[rxn].val.err / dGr_dict[rxn].val.val > variation_tres || abs(dGr_dict[rxn].val.val) == Inf
                    mod_bounds[rxn] = (0.0, 0.0)
                end
            catch
                mod_bounds[rxn] = (0.0, 0.0)
            end
        end
    end

    return database, mod_bounds
end
