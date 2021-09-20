
"""
    function collect_dGf(database; ph::Float64 = 7.0, pmg::Float64 = 2.0, i_strengh::Quantity = 100.0u"mM", temp::Quantity = 25u"°C", db_id = bigg)
Collects ΔG of formation for all metabolites in `database` via Component Contribution.
Kwargs setup eQuilibrator conditions.
"""
function collect_dGf(
    database;
    ph::Float64 = 7.0,
    pmg::Float64 = 2.0,
    i_strengh::Quantity = 100.0u"mM",
    temp::Quantity = 25u"°C",
    db_id = bigg,
)
    println("eQuilibrator is initialized. Please wait...")
    equil = eQuilibrator.Equilibrator(
        pH = ph,
        pMg = pmg,
        ionic_strength = i_strengh,
        temperature = temp,
    )
    dGf_dict = OrderedDict()
    n = length(metabolites(database))
    p = Progress(n, 1)
    for met in metabolites(database)
        rxn_str = string(" = ", met[1:end-2])
        try
            @suppress begin
                dGf_dict[met] =
                    physiological_dg_prime(equil, db_id(rxn_str); balance_warn = false)
            end
        catch
            dGf_dict[met] = "NA"
        end
        next!(p)
    end
    return dGf_dict
end

"""
    function collect_dGr_bounds(database; return_opts::Int = 3, ph::Float64 = 7.0, pmg::Float64 = 2.0, i_strengh::Quantity = 100.0u"mM", temp::Quantity = 25u"°C", db_id = bigg)
Collects ΔG of reaction for every reaction in `database` via Component Contribution.
Returns dictionary of ΔG of reactions (`return_opts = 1`), the calculated bounds for MinFlux/MinRxn (`return_opts = 2`) or both (`return_opts = 3`).
Kwargs setup eQuilibrator conditions and return.
"""
function collect_dGr_bounds(
    database;
    return_opts::Int = 3,
    ph::Float64 = 7.0,
    pmg::Float64 = 2.0,
    i_strength::Quantity = 100.0u"mM",
    temp::Quantity = 25u"°C",
    db_id = bigg,
)
    println("eQuilibrator is initialized. Please wait...")
    equil = eQuilibrator.Equilibrator(
        pH = ph,
        pMg = pmg,
        ionic_strength = i_strength,
        temperature = temp,
    )

    R = 8.31446261815324u"J/(K*mol)"
    R = uconvert(u"kJ/(K*mol)", R)


    println("Calculating standard ΔG of reaction:")
    dGr_standard_dict = OrderedDict()
    n = length(reactions(database))
    p = Progress(n, 1)
    for rxn in reactions(database)
        rxn_string = stoichiometry_string(
            reaction_stoichiometry(database, rxn),
            format_id = x -> x[1:end-2],
        )
        try
            @suppress begin
                dGr_standard_dict[rxn] =
                    standard_dg_prime(equil, db_id(rxn_string), balance_warn = false)
            end
        catch
            dGr_standard_dict[rxn] = "NA"
        end
        next!(p)
    end


    println("Calculating flux bounds:")
    dGr_bounds = OrderedDict()
    n = length(reactions(database))
    p2 = Progress(n, 1)

    for (rxn, dG) in dGr_standard_dict
        try
            @suppress begin
                reactants = [v for (k, v) in reaction_stoichiometry(database, rxn) if v > 0]
                products =
                    [abs(v) for (k, v) in reaction_stoichiometry(database, rxn) if v < 0]
                qmax =
                    prod(fill(0.1, length(products)) .^ products) /
                    prod(fill(1 * 10^-6, length(reactants)) .^ reactants)
                qmin =
                    prod(fill(1 * 10^6, length(products)) .^ products) /
                    prod(fill(0.1, length(reactants)) .^ reactants)
                u_calc = R * uconvert(u"K", temp) * log(qmax)
                l_calc = R * uconvert(u"K", temp) * log(qmin)
                println(typeof(dG))
                dGu = dG.val + u_calc.val
                dGl = dG.val + l_calc.val
                dGr_bounds[rxn] = (dGl, dGu)
            end
        catch
            dGr_bounds[rxn] = "NA"
        end
        next!(p2)
    end
    if return_opts == 3
        return dGr_standard_dict, dGr_bounds
    elseif return_opts == 2
        return dGr_bounds
    elseif return_opts == 1
        return dGr_standard_dict
    else
        println("Return options not specified. Switched to complete return.")
        return dGr_standard_dict, dGr_bounds
    end
end


function reaction_bounds(bound_dc; M = 1000)
    binary_dc_lb = OrderedDict()
    binary_dc_ub = OrderedDict()
    for (rxn, vals) in bound_dc
        #LB binaries
        if vals[1] < 0
            LB = 0
        else
            LB = -M
        end
        if vals[2] < 0
            UB = M
        else
            UB = 0
        end
        binary_dc_lb[rxn] = LB
        binary_dc_ub[rxn] = UB
    end
    return binary_dc_lb, binary_dc_ub
end



