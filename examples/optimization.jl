"""
Test script forOptStoic/MinFlux procedure. 
Model construction and package imports can be found in `testmodel_construction.jl`.
Output is expected to be as follows:
Dict{Any, Any} with 10 entries:
  "r5"      => -1.0
  "r2"      => 1.0
  "EX_m6"   => 2.0
  "biomass" => -2.0
  "EX_co1"  => -2.0
  "EX_m1"   => -1.0
  "EX_h"    => 2.0
  "EX_co2"  => 2.0
  "r3"      => 1.0
  "r4"      => 1.0
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
    m, mf = build_MinFlux_model(m, os_solution, minflux_bounds, optimizer = Gurobi)
    mf = MinFlux(mf)
    MinFlux_solution(m, mf)
end
