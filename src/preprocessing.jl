"""
    reaction_dg_bounds(database, dgr_dict, low_con::Float64, high_con::Float64)

Calculates the minimal and maximal transformed ΔG of reaction for every ΔG'° depending on assumed low and high concentration values (`low_con` and `high_con`, respectively). Returns a dictionary with tuples of those values.
The standard transformed ΔGs of reaction are passed via `dgr_dict`.
"""
function reaction_dg_bounds(
        database, 
        dgr_dict, 
        low_concentration::Float64, 
        high_concentration::Float64; 
        temperature::Float64 = 298.15,
    )
    
    ln_low = log(low_concentration)
    ln_high = log(high_concentration)

    R = 8.31446261815324u"J/(K*mol)"
    R = uconvert(u"kJ/(K*mol)", R)

    dg_bounds = OrderedDict()
    for (rxn, dg) in dgr_dict
        educts = [abs(v) for (k, v) in reaction_stoichiometry(database, rxn) if v < 0]
        products = [abs(v) for (k, v) in reaction_stoichiometry(database, rxn) if v > 0]
        qmin = ln_low * sum(products) - ln_high * sum(educts)
        qmax = ln_high * sum(products) - ln_low * sum(educts)
        dg_min = dg.val.val + R.val * temperature * qmin # maximum value
        dg_max = dg.val.val + R.val * temperature * qmax # minimum value
        
        dg_bounds[rxn] = (dg_min, dg_max)
    end
    return dg_bounds
end

"""
    flux_bounds(bound_dc; M = 1000)
        
Returns lower and upper flux bounds for MinFlux procedure. M is set to be a high number (default at 1000). 
The energy bounds can be calculated via the reaction_dg_bounds function.
"""
function flux_bounds(energy_bounds; M = 1000)
    binary_dc = OrderedDict()
    for (rxn, vals) in energy_bounds
        # Calculating binaries
        if vals[2] ≤ 0
            LB = 0
        else
            LB = -M
        end
        if vals[1] < 0
            UB = M
        else
            UB = 0
        end
        if LB > UB
            LB = 0
            UB = 0
        end
        binary_dc[rxn] = (LB, UB)
    end
    return binary_dc
end

"""
    MinFlux_preprocessing(
        model,
        standard_Gibbs_dictionary;
        lower_concentration = 1.0e-6,
        upper_concentration = 1.0e-1,
        temperature = 298.15,
        M = 1000,
    )

Complete preprocessing routine for constraining data for the MinFlux procedure. Returns a dictionary with lower and upper bounds for the corresponding optimization constraints. Keyword arguments describe the default values for the lower and upper concentrations (lower_concentration, upper_concentration) for the thermodynamic data, the system temperature (temperature) and the flux margin (M).
"""
function MinFlux_preprocessing(
        model,
        standard_Gibbs_dictionary;
        lower_concentration = 1.0e-6,
        upper_concentration = 1.0e-1,
        temperature = 298.15,
        M = 1000,
    )
    dgs = reaction_dg_bounds(model, standard_Gibbs_dictionary, lower_concentration, upper_concentration, temperature = temperature)
    return flux_bounds(dgs, M = M)
end
