"""
Thermdynamic data collection

Function names list (order by appearance in code):
    metacyc
    metacycIllParsed
    collect_formation_dg (internal eQuilibrator setup)
    collect_formation_dg (overload without internal eQuilibrator setup)
    collect_reaction_dg
    reaction_dg_bounds
    collect_dGr_bounds
    reaction_bounds
"""

"""
    $(SIGNATURES)
Returns reaction string with a prefixed "metacyc.compound:".
"""
function metacyc(str)
    eQuilibrator._parse_reaction_string(str, "metacyc.compound")
end

"""
    $(SIGNATURES)
Returns modified string with prefixed "metacyc.compound:". If the MetaCyc-based model was created via the Python-based "Moped" module, metabolite and reaction names may have artifacts ("__45__" and "__45__" substrings) that should be replaced by a hyphen or periodt, respectively, before ΔG estimations.
"""
function metacycIllParsed(str)
    nstr = replace(str, "__45__46__45__" => "-46-")
    nstr = replace(nstr, "__46__" => ".")
    nstr = replace(nstr, "__45__" => "-")
    nstr = replace(nstr, "__43__" => "+")
    eQuilibrator._parse_reaction_string(nstr, "metacyc.compound")
end

"""
    collect_formation_dg(
        database,
        db_id;
        dGf_dict = OrderedDict(),
        ph::Float64 = 7.0,
        pmg::Float64 = 2.0,
        i_strengh::Quantity = 100.0u"mM",
        temp::Quantity = 25u"°C",
    )

Collects ΔG of formation for all metabolites in `database` via Component Contribution. The data bank source must be passed by `db_id` (Options: `bigg`, `kegg`, `metanex`, `chebi`). Alternatively a dictionary of ΔGs of formation `dGf_dict` can be passed to extend it by new entries for the database metabolites.
Keyword arguments setup eQuilibrator conditions.

# Example
```
collect_formation_dg(
    ecoli_core,
    bigg
)
```
"""
function collect_formation_dg(
    database,
    db_id;
    dGf_dict = OrderedDict(),
    met_prefix = "",
    ph::Float64 = 7.0,
    pmg::Float64 = 2.0,
    i_strengh::Quantity = 100.0u"mM",
    temp::Quantity = 25u"°C",
)
    println("eQuilibrator is initialized. Please wait...")

    equil = eQuilibrator.Equilibrator(
        pH = ph,
        pMg = pmg,
        ionic_strength = i_strengh,
        temperature = temp,
    )

    n = length(metabolites(database))
    p = Progress(n, 1)
    for met in metabolites(database)
        if met ∉ keys(dGf_dict)
            rxn_str = string(" = ", met[length(met_prefix)+1:end-2])
            try
                @suppress begin
                    dGf_dict[met] =
                        standard_dg_prime(equil, db_id(rxn_str); balance_warn = false)
                end
            catch
                dGf_dict[met] = "NA"
            end
        end
        next!(p)
    end
    return dGf_dict
end

"""
    collect_formation_dg(
        database,
        equil,
        db_id;
        met_prefix = "",
        dGf_dict = OrderedDict(),
    )

Collects ΔG of formation for all metabolites in `database` via Component Contribution. The data bank source must be passed by `db_id` (Options: `bigg`, `kegg`, `metanex`, `chebi`). An eQuilibrator setup has to be passed by `equil` (see also: https://github.com/stelmo/Equilibrator.jl). Alternatively a dictionary of ΔGs of formation `dGf_dict` can be passed to extend it by new entries for the database metabolites.

# Example
```
collect_formation_dg(
    ecoli_core,
    eQuilibrator.Equilibrator(ionic_strength=150.0u"mM"),
    bigg
)
```
"""
function collect_formation_dg(
    database,
    equil,
    db_id;
    met_prefix = "",
    dGf_dict = OrderedDict(),
)
    println("ΔG of formation are collected.")
    n = length(metabolites(database))
    p = Progress(n, 1)
    for met in metabolites(database)
        if met ∉ keys(dGf_dict)
            rxn_str = string(" = ", met[length(met_prefix)+1:end-2])
            try
                @suppress begin
                    dGf_dict[met] =
                        standard_dg_prime(equil, db_id(rxn_str); balance_warn = false)
                end
            catch
                dGf_dict[met] = measurement(1000000u"kJ/mol", 10000000u"kJ/mol") #setting infeasible value
            end
        end
        next!(p)
    end
    return dGf_dict
end

"""
function collect_reaction_dg(
    database,
    equil,
    db_id;
    met_prefix = "",
    dgr_dict = OrderedDict(),
)
Collects standard transformed ΔG of reactions for every reaction in an SBML model. The eQuilibrator setup (see `eQuilibrator.Equilibrator()`) is passed as the `equil` argument. The source database (e.g. BiGG) is passed by `db_id`. If a dictionary with ΔGs of reaction already exist and should be expand, pass the dictionary to the `dgr_dict` keyword argument.
"""
function collect_reaction_dg(
    database,
    equil,
    db_id;
    met_prefix = "",
    dgr_dict = OrderedDict(),
)
    println("Standard ΔG of reactions are collected.")
    n = length(metabolites(database))
    p = Progress(n, 1)
    for rxn in reactions(database)
        if rxn ∉ keys(dgr_dict)
            rxn_string = stoichiometry_string(
                reaction_stoichiometry(database, rxn),
                format_id = x -> db_id(x[length(met_prefix)+1:end-2]),
            )
            try
                @suppress begin
                    dgr_dict[rxn] = standard_dg_prime(equil, rxn_string, balance_warn = false)
                end
            catch
                dgr_dict[rxn] = measurement(1000000u"kJ/mol", 10000000u"kJ/mol") #setting infeasible value
            end
        end
        next!(p)
    end
    return dgr_dict
end

"""
    collect_dGr_bounds(
        database;
        return_opts::Int = 3,
        ph::Float64 = 7.0,
        pmg::Float64 = 2.0,
        i_strength::Quantity = 100.0u"mM",
        temp::Quantity = 25u"°C",
        db_id = bigg,
        met_suffix = "",
    )

DEV NOTE: Necessary for MinFlux/MinRxn procedure, currently work in progress. Function might change severely.
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
    met_prefix = "",
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
            format_id = x -> x[length(met_prefix)+1:end-2],
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
                reactants = 
                    [abs(v) for (k, v) in reaction_stoichiometry(database, rxn) if v < 0]
                products =
                    [abs(v) for (k, v) in reaction_stoichiometry(database, rxn) if v > 0]
                qmax =
                    prod(fill(0.1, length(products)) .^ products) /
                    prod(fill(1 * 10^-6, length(reactants)) .^ reactants)
                qmin =
                    prod(fill(1 * 10^-6, length(products)) .^ products) /
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
    if return_opts == 3 # Returns both ΔG of reaction dictionary and the calculated bounds
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

function reduce_model(database, dg_dictionary)
    for rxn in reactions(database)
        if rxn ∉ keys(dg_dictionary)
            database = remove_reaction(database, rxn)
        end
    end
    return database
end