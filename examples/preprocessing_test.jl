include("testmodel_construction.jl")
flush(stdout)
m # model
reaction_dG # ΔG'° of reaction dictionary
minflux_bounds = MinFlux_preprocessing(m, reaction_dG)