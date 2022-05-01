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
    substrates::Vector{String},
    targets::Vector{String},
    co_reactants::Vector{String},
    optimizer;
    db_id = bigg,
	energy_dc = Dict(),
    variable_bounds = Dict("substrates" => -1, "reactants" => (-15, 15)),
    custom_bounds = Dict(),
    dG_thres = -5,
    IP = true,
) where {SM<:COBREXA.StandardModel}

    # Constructing the LP model
    os_model = Model(optimizer.Optimizer)

    # Adds custom bounds to the variable bounds
    merge!(variable_bounds, custom_bounds)

    # Collecting elements for constraints
    participants = vcat(substrates, targets, co_reactants)
    elements = Vector{String}()
    for met in participants
        elements = vcat(elements, collect(keys(metabolite_formula(database, met))))
    end
    elements = Set(elements)

    # Metabolite formulae are extended by "empty" elements to allow mass balance in constraint construction.
    ### INEFFICIENT. WILL BE CHANGED SOON ###
    ext_formulae = _extend_formulae(database)

    # Assigning ΔGs for constraints
    dgf_substrates = Vector()
    for met in substrates
        push!(dgf_substrates, energy_dc[met].val.val)
        if met ∉ keys(variable_bounds)
            variable_bounds[met] = variable_bounds["substrates"]
        end
    end

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
            S[it = 1:length(substrates)] == variable_bounds[substrates[it]],
            (start = it, base_name = substrates[it], integer = IP)
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
        sum(S[it] * dgf_substrates[it] for it = 1:length(substrates)) 
        ≤ dG_thres
    )
    @constraint(
        os_model,
        ChargeBalance,
        sum(T[it] * metabolite_charge(database, targets[it]) for it = 1:length(targets)) +
        sum(
            CoR[it] * metabolite_charge(database, co_reactants[it]) for
            it = 1:length(co_reactants)
        ) +
        sum(S[it] * metabolite_charge(database, substrates[it]) for it = 1:length(substrates))
        == 0
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
            sum(S[it] * ext_formulae[substrates[it]][elem] for it = 1:length(substrates))
            == 0
        )

        c_norm = -1 * sum(variable_bounds[substrates[it]] * ext_formulae[substrates[it]]["C"] for it = 1:length(substrates)) # carbon normalization factor
        @objective(os_model, Max, sum(T) / c_norm) # Objective function
    end

    return os_model
end

"""
    build_MinFlux_model(
        database,
        optstoic_solution,
        dGr_dict = Dict();
        L = [],
        optimizer,
        loopless = false,
        N_red = [],
        ϵ = 10^(-5),
        M = 1000
    )
MinFlux formulation. Returns MinFlux MILP model. Takes an overall stoichiometry or a previously calculated OptStoic solution dictionary (`optstoic_solution`) as argument. Stoichiometry may be structured as a dictionary with metabolite IDs as keys and the corresponding stoichiometric coefficients as values (negative for substrates, positive for products). The IDs need to match the corresponding metabolite IDs in the metabolic model (`database`). The flux_bounds_dictionary argument is a dictionary with the preprocessed boundaries (see function `flux_bounds` and `MinFlux_preprocessing`). Loopless constraints can be added by setting `loopless` as `true`. A precomiled matrix `N_red` for loopless constraints can additionally be passed to avoid automatic matrix calculation. 
)
"""
function build_MinFlux_model(
    database,
    optstoic_solution,
    flux_bounds_dictionary;
    L = [],
    optimizer,
    loopless = false,
    N_red = [],
    ϵ = 10^(-5),
    M = 1000
)
    mf_model = Model(optimizer.Optimizer)

    # If loopless constraints are activated and no constraint matrix L is passed via keyword, a matrix is constructed locally.
    if loopless
        if isempty(N_red)
            N_red = Nred(database)
        end
        red_dimensions = size(N_red)[1]
    end

    ex_vec = Vector()
    for (comp, coeff) in optstoic_solution
        push!(ex_vec, Reaction(string("EX_", comp), Dict(comp => -1)))
        flux_bounds_dictionary[string("EX_", comp)] = (coeff, coeff)
    end
    for ex in ex_vec
        if ex ∉ reactions(database)
            add_reaction!(database, ex)
        end
    end

    lbs = Vector()
    ubs = Vector()
    for x in collect(values(flux_bounds_dictionary))
        push!(lbs, x[1]) # Collecting lower bounds
        push!(ubs, x[2]) # Collecting upper bounds
    end

    rxn_size = length(reactions(database))
    @variable(mf_model, V[i = 1:rxn_size], Int)
    @variable(mf_model, X[i = 1:rxn_size])
    @variable(mf_model, Vf[i = 1:rxn_size])
    @variable(mf_model, Vr[i = 1:rxn_size])
    @variable(mf_model, Yf[i = 1:rxn_size], Bin)
    @variable(mf_model, Yr[i = 1:rxn_size], Bin)

    # Adding loopless constraints, if loopless is true
    if loopless
        @variable(mf_model, G[i = 1:red_dimensions])
        @variable(mf_model, a[i = 1:red_dimensions], Bin)

        @constraint(mf_model, LoopBalance, N_red' * G .== 0)
        @constraint(mf_model, L1, G .>= -M * a + (1 .- a))
        @constraint(mf_model, L2, G .>= -a .+ M * (1 .- a))
        @constraint(mf_model, L3, [V[it] for it in 1:length(V) if !startswith(reactions(database)[it], "EX_")] .>= -M * (1 .- a))
        @constraint(mf_model, L4, [V[it] for it in 1:length(V) if !startswith(reactions(database)[it], "EX_")] .<= M * a)
    end

    @constraint(mf_model, MassBalance, stoichiometry(database) * V .== balance(database))
    @constraint(mf_model, LowerBounds, lbs .<= V)
    @constraint(mf_model, Upperbounds, ubs .>= V)
    @constraint(mf_model, FluxDiff, V .== Vr - Vf)
    @constraint(mf_model, FluxSum, X .== Vr + Vf)
    @constraint(mf_model, ForwardFlux, Vf .>= 0)
    @constraint(mf_model, ReverseFlux, Vr .>= 0)
    @constraint(mf_model, abscon1, X .>= V)
    @constraint(mf_model, abscon2, X .>= -V)
    @constraint(mf_model, Vf .>= ϵ * Yf)
    @constraint(mf_model, M * Yf .>= Vf)
    @constraint(mf_model, Vr .>= ϵ * Yr)
    @constraint(mf_model, M * Yf .>= Vr)
    @constraint(mf_model, BinActive, 1 .>= Yf + Yr)
    if !isempty(L)
        for k in 1:length(L)
            @constraint(mf_model, base_name = string("IntegerCut_it", k),
            sum([1-Yf[it]-Yr[it] for it in 1:rxn_size if sum(L[k][it]) == 1 && !startswith(reactions(database)[it], "EX_")]) >= 1)
        end
    end
    @objective(mf_model, Min, sum(X))

    return database, mf_model

end


"""
    extract_binaries(minflux_model, arr = [], ϵ = 0)

Creates tuple vector for cut constraints. Can be added to previous vector (`arr`) to extend previous cut constraints. `ϵ` is used as a threshold for residual float values in the binary LP variables.
"""
function extract_binaries(minflux_model, arr = [], ϵ = 0)
    new_path = []
    for i in 1:length(minflux_model[:Yf])
        if value(minflux_model[:Vf][i]) <= 2*ϵ
            yf = 0.0
        else
            yf = value(minflux_model[:Yf][i])
        end
        if value(minflux_model[:Vr][i]) <= 2*ϵ
            yr = 0.0
        else
            yr = value(minflux_model[:Yr][i])
        end

        push!(new_path, (yr, yf))
    end

    return push!(arr, new_path)
end