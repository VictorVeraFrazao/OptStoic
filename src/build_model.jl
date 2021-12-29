"""
LP model construction routines

Function names list (order by appearance in code):
    build_OptStoic_model
    build_MinFlux_model
    build_MinRxn_model
"""

"""
    build_OptStoic_model(
        database::S,
        substrate::String,
        targets::Vector{String},
        co_reactants = String[],
        energy_dc = Dict(),
        optimizer;
        energy_dc = Dict(),
        db_id = bigg,
        variable_bounds = Dict("substrate" => -1, "reactants" => (-15, 15)),
        custom_bounds = Dict(),
        dG_thres = -5,
        IP = true,
    ) where {S<:COBREXA.StandardModel}

Construction of LP model for OptStoic procedure. OptStoic procedure is formulated via this optimisation problem:
```
max f(scⱼ)
s.t.
    ∑(nᵤⱼ scⱼ) = 0, ∀ u∈U   (Mass Balance)
    ∑(eⱼ scⱼ) = 0,          (Charge Balance)
    ∑(ΔGⱼ scⱼ) ≤ ΔGᵐⁱⁿ      (Energy Balance)
```
This formulation was derived and shortened for display. See "Chowdhury, A., Maranas, C. Designing overall stoichiometric conversions and intervening metabolic reactions. Sci Rep 5, 16009 (2015). https://doi.org/10.1038/srep16009" for more information. It is possible to formulate this problem as IP or MIP, hence the keyword binary `IP` with `true` for IP and `false` for MIP.

The `substrate` should be the limiting carbon source of the desired biochemical conversion (e.g. glucose). The `targets` include every target metabolite that should be maximised. The `co_reactants` include every metabolite that is allowed to appear in the stoichiometry but is not supposed to be targeted for maximisation. The `optimizer` has to be set to a `JuMP`-compatible optimiser, such as `JuMP.Gurobi` or `JuMP.GLPK` and supports (MI)LP (see https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers for compatibility issues).

The construction process includes the collecting process of ΔG of formation for the metabolites within the `database`. Optionally a dictionary with the necessary thermodynamic data can be passed to the function via `energy_dc`. The default limits for the stoichiometric coefficients are set in `variable_bounds` and can be altered. Custom bounds for specific metabolites can be passed by the `custom_bounds` dictionary. The ΔGᵐⁱⁿ (`dG_thres`) is set to -5 kJ/mol by default.
# Example
```
Glucose -> Acetate conversion (anaerobic, E. coli core model (from http://bigg.ucsd.edu/models/e_coli_core))
build_OptStoic_model(
    ecoli_core,
    "glc__D_e",
    ["ac_e"],
    ["co2_e", "o2_e", "h2o_e", "h_e"],
    Gurobi,
    custom_bounds = Dict("o2_e" => (0, 0))
)
```
"""
function build_OptStoic_model(
    database::SM,
    substrate::String,
    targets::Vector{String},
    co_reactants::Vector{String},
    optimizer;
    db_id = bigg,
	energy_dc = Dict(),
    variable_bounds = Dict("substrate" => -1, "reactants" => (-15, 15)),
    custom_bounds = Dict(),
    dG_thres = -5,
    IP = true,
) where {SM<:COBREXA.StandardModel}

    # Constructing the LP model
    os_model = Model(optimizer.Optimizer)

    # Collecting ΔG of formation if no data is passed manually. Throws error if the given data does not cover the metabolite data of the given model
    if isempty(energy_dc)
        energy_dc = collect_formation_dg(database, db_id)
    else
        for met in metabolites(database)
            if met ∉ keys(energy_dc)
                throw(MissingValueError(met, "metabolite ΔG or formation not defined"))
            end
        end
    end

    merge!(variable_bounds, custom_bounds) # Adds custom bounds to the variable bounds

    # Collecting elements for constraints
    participants = vcat(targets, co_reactants)
    push!(participants, substrate)
    elements = Vector{String}()
    for met in participants
        elements = vcat(elements, collect(keys(metabolite_formula(database, met))))
    end
    elements = Set(elements)

    # Metabolite formulae are extended by empty values to allow mass balance in constraint construction.
    # DEV NOTE: Process may be inefficient. Implementation of alternate processes recommended.
    ext_formulae = _extend_formulae(database)

    # Assigning ΔGs for constraints
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
            S == variable_bounds["substrate"], (base_name = substrate, integer = IP)
            1 ≤ T[it = 1:length(targets)] ≤ variable_bounds[targets[it]][2],
            (start = it, base_name = targets[it], integer = IP)
            variable_bounds[co_reactants[it]][1] ≤
            CoR[it = 1:length(co_reactants)] ≤
            variable_bounds[co_reactants[it]][2],
            (start = it, base_name = co_reactants[it], integer = IP)
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

        @objective(os_model, Max, sum(T) / -1 * S) # Objective function
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
MinFlux optimisation model construction. WORK IN PROGRESS!!! Documentation will follow.
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

    return database, mf_model

end

"""
    function build_MinRxn_model(
        database,
        optstoic_solution,
        dGr_dict = Dict();
        optimizer,
    )
MinRxn optimisation model construction. WORK IN PROGRESS!!! Documentation will follow.
"""
function build_MinRxn_model(
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
    @variable(mf_model, Y[i = 1:length(reactions(database))], Bin)
    @constraint(mf_model, MassBalance, stoichiometry(database) * V .== balance(database))
    @constraint(mf_model, LowerBounds, lbs .* Y .<= V)
    @constraint(mf_model, Upperbounds, ubs .* Y .>= V)
    for it in 1:length(Y)
        @constraint(mf_model, Y[it] => {lbs[it] <= V[it] <= ubs[it]})
        @constraint(mf_model, !Y[it] => {V[it] == 0})
    end
    @objective(mf_model, Min, sum(Y))

    return mf_model

end