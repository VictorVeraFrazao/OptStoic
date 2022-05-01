"""
    Quick MinFlux preprocessing test for MinFlux constraints.
"""

include("testmodel_construction.jl")
m # model
reaction_dG # ΔG'° of reaction dictionary
minflux_bounds = MinFlux_preprocessing(m, reaction_dG)
