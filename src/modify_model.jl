"""
    extend_formulae(database)
Adapts stoichiometric formulae of a `MetabolicModel` for OptStoic constraints.
Returns a dictionary of updated formulae of the model, so that every present metabolite has an entry for all chemical elements within the models.
"""
function extend_formulae(database)
    updated_formulae = Dict()
    elem_set = Vector()
    for m in metabolites(database)
        elem_set = vcat(elem_set, collect(keys(metabolite_formula(database, m))))
    end
    elem_set = Set(elem_set)
    for m in metabolites(database)
        formula = metabolite_formula(database, m)
        t = setdiff(elem_set, Set(collect(keys(formula))))
        for i in t
            if i ∉ keys(formula)
                formula[i] = 0
            end
        end
        updated_formulae[m] = formula
    end
    return updated_formulae
end

"""
    function adjust_model(database, variation_tres =.5, dGr_dict =Dict())
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
