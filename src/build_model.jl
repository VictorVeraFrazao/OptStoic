"""
    build_OptStoic_model(
        database,
        substrate::String,
        targets::Vector{String},
        co_reactants = String[],
        energy_dc = Dict();
        optimizer;
        variable_bounds = Dict("substrate" => 1, "reactants" => 15),
        dG_thres = -5,
    )

Returns model for OptStoic procedure. Database should be of the `COBREXA.StandardModel` type. The `substrate` string should be the ID of the limiting metabolite (e.g. glucose). The `targets` vector should only contain the products that should be contained. If metabolites should be allowed as co-factors or co-reactants in the overall stoichiometry, pass them via the `co_reactants` vector. The ΔG of formation for the database metabolites can be passed via `energy_dc` as a dictionary. If none is passed they are calculated automatically (be aware that this increases runtime). Kwargs are used for constraint bounds.
    
    Work in progress, solver needs to be able to handle MIPs? exclusively.

# Example
```
build_OptStoic_model(
    ecoli_core,
    "glc__D_e",
    ["pyr_e", "atp_c"],
    ["adp_c", "nadh_c", "nad_c", "h2o_e", "h_e"],
    optimizer = Gurobi
)
```
"""
function build_OptStoic_model(
    database,
    substrate::String,
    targets::Vector{String},
    co_reactants::Vector{String},
    optimizer;
	energy_dc = Dict(),
    variable_bounds = Dict("substrate" => 1, "reactants" => 15),
    dG_thres = -5,
)
    os_model = Model(optimizer.Optimizer)

    if isempty(energy_dc)
        energy_dc = collect_dGf(database)
    end

    # Collecting elements for constraints
    participants = vcat(targets, co_reactants)
    push!(participants, substrate)
    elements = Vector{String}()
    for met in participants
        elements = vcat(elements, collect(keys(metabolite_formula(database, met))))
    end
    elements = Set(elements)

    ext_formulae = extend_formulae(database)
    dgf_substrate = energy_dc[substrate].val.val

    dgf_targets = Vector()
    for met in targets
        push!(dgf_targets, energy_dc[met].val.val)
        if met ∉ keys(variable_bounds)
            variable_bounds[met] = variable_bounds["reactants"]
        end
    end

    dgf_cor = Vector()
    for met in co_reactants
        push!(dgf_cor, energy_dc[met].val.val)
        if met ∉ keys(variable_bounds)
            variable_bounds[met] = variable_bounds["reactants"]
        end
    end

    # Assigning variables and constraints
    @variables(
        os_model,
        begin
            S == -1 * variable_bounds["substrate"], (base_name = substrate, integer = true)
            1 ≤ T[it = 1:length(targets)] ≤ variable_bounds[targets[it]],
            (start = it, base_name = targets[it], integer = true)
            -1 * variable_bounds[co_reactants[it]] ≤
            CoR[it = 1:length(co_reactants)] ≤
            variable_bounds[co_reactants[it]],
            (start = it, base_name = co_reactants[it], integer = true)
        end
    )

    @constraint(
        os_model,
        EnergyBalance,
        sum(T[it] * dgf_targets[it] for it = 1:length(targets)) +
        sum(CoR[it] * dgf_cor[it] for it = 1:length(co_reactants)) +
        S * dgf_substrate ≤ dG_thres
    )
    @constraint(
        os_model,
        ChargeBalance,
        sum(T[it] * metabolite_charge(database, targets[it]) for it = 1:length(targets)) +
        sum(
            CoR[it] * metabolite_charge(database, co_reactants[it]) for
            it = 1:length(co_reactants)
        ) +
        S * metabolite_charge(database, substrate) == 0
    )
    for elem in elements
        @constraint(
            os_model,
            (base_name = string("MassBalance_", elem)),
            sum(T[it] * ext_formulae[targets[it]][elem] for it = 1:length(targets)) +
            sum(
                CoR[it] * ext_formulae[co_reactants[it]][elem] for
                it = 1:length(co_reactants)
            ) +
            S * ext_formulae[substrate][elem] == 0
        )

        @objective(os_model, Max, sum(T) / -1 * S)
    end

    return os_model
end

"""
    function build_MinFlux_model(
        database,
        optstoic_solution,
        dGr_dict;
        optimizer,
    )
Builds model for MinFlux procedure. The given `database` is automatically adapted and the ΔG of reactions are calculated internally. Alternatively a manually adapted model, together with a ΔG list/dictionary can be hand over (WORK IN PROGRESS!).
"""
function build_MinFlux_model(
    database,
    optstoic_solution,
    dGr_dict = Dict(); 
    optimizer,
)

    if isempty(dGr_dict)
        database, dGr_dict = adjust_model(database)
    end
    b_dict = reaction_bounds(dGr_dict)
    mf_model = Model(optimizer.Optimizer)

    ex_vec = Vector()

    for (comp, coeff) in optstoic_solution
        push!(ex_vec, Reaction(string("EX_", comp), Dict(comp => -1)))
        b_dict[string("EX_", comp)] = (coeff, coeff)
    end
    for ex in ex_vec
        add_reaction!(database, ex)
    end

    lbs = Vector()
    ubs = Vector()
    for x in collect(values(b_dict))
        push!(lbs, x[1]) # Collecting lower bounds
        push!(ubs, x[2]) # Collecting upper bounds
    end

    @variable(mf_model, V[i = 1:length(reactions(database))])
    @variable(mf_model, X[i = 1:length(reactions(database))])
    @constraint(mf_model, MassBalance, stoichiometry(database) * V .== balance(database))
    @constraint(mf_model, LowerBounds, lbs .<= V)
    @constraint(mf_model, Upperbounds, ubs .>= V)
    @constraint(mf_model, abscon1, X .>= V)
    @constraint(mf_model, abscon2, X .>= -V)
    @objective(mf_model, Min, sum(X))

    return mf_model

end
