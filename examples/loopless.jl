"""
Test script for loopless constraints in MinFlux. 
Model construction and package imports can be found in `testmodel_construction.jl`.
"""

include("testmodel_construction.jl")

minflux_bounds = MinFlux_preprocessing(m, reaction_dG)
substrates = ["m1"]
targets = ["m6"]
cofactors = ["co1", "co2", "h"] #allow cofactors

bounds = Dict("substrates" => -1, "reactants" => (-15, 15))

os = build_OptStoic_model(m, substrates, targets, cofactors, Gurobi, energy_dc = formation_dG, variable_bounds = bounds)
os = optimize_OptStoic(os)
os_solution = solution_OptStoic(os)

if !isnothing(os_solution)
    m, mf = build_MinFlux_model(m, os_solution, minflux_bounds, optimizer = Gurobi, loopless = true)
    mf = MinFlux(mf)
    MinFlux_solution(m, mf)
end #Output should be the same as for `optimization.jl`