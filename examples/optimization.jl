using Pkg
Pkg.activate(joinpath(Pkg.devdir(), "OptStoic/"))
using OptStoic
using Gurobi
using JuMP
using COBREXA
using CSV
using Measurements, Unitful

"""
Construction of model
"""
m = StandardModel("SmallModel")
m1 = Metabolite("m1", charge = 0, formula = "C6H12O6")
m2 = Metabolite("m2", charge = -2, formula = "C6H11O9P1")
m3 = Metabolite("m3", charge = -2, formula = "C6H11O9P1")
m4 = Metabolite("m4", charge = -4, formula = "C6H10O12P2")
m5 = Metabolite("m5", charge = -2, formula = "C3H5O6P1")
m6 = Metabolite("m6", charge = -2, formula = "C3H5O6P1")
co1 = Metabolite("co1", charge = -3, formula = "C10H12N5O13P3")
co2 = Metabolite("co2", charge = -2, formula = "C10H12N5O10P2")
h = Metabolite("h", charge = 1, formula = "H1")

@add_reactions! m begin
    "r1", nothing → m1, 0, 100
    "r2", m1 + co1 → m2 + co2 + h, 0, 100
    "r3", m2 → m3, 0, 100
    "r4", m3 + co1 → m4 + co2 + h, 0, 100
    "r5", m5 + m6 → m4, 0, 100
    "r6", m5 → m6, -100, 100
    "r7", nothing → co1, 0, 100
    "r8", co2 → nothing, 0, 100
    "r9", nothing → h, -100, 100
    "biomass", m6 → nothing, 0, 100
end

gs = [Gene("g$i") for i = 1:5]

m.reactions["biomass"].objective_coefficient = 1.0

add_genes!(m, gs)
add_metabolites!(m, [m1, m2, m3, m4, m5, m6, co1, co2, h])


formation_dG = Dict(
# Dictionary with ΔG of formation for every metabolite
    "m1" => measurement(-404.35, 1.13)u"kJ/mol",
    "m2" => measurement(-1293.33, 1.89)u"kJ/mol",
    "m3" => measurement(-1293.15, 1.27)u"kJ/mol",
    "m4" => measurement(-2187.92, 2.02)u"kJ/mol",
    "m5" => measurement(-1085.78, 0.56)u"kJ/mol",
    "m6" => measurement(-1080.17, 0.67)u"kJ/mol",
    "co1" => measurement(-2277.33, 1.49)u"kJ/mol",
    "co2" => measurement(-1402.53, 1.23)u"kJ/mol",
    "h" => measurement(0, 0)u"kJ/mol",
)

reaction_dG = Dict(
    # Dictionary with ΔG of reactions for every reaction
    "r1" => measurement(-404.35, 1.13)u"kJ/mol",
    "r2" => measurement(-14.177, 1.92)u"kJ/mol",
    "r3" => measurement(-1.14, 1.63)u"kJ/mol",
    "r4" => measurement(-19.95, 1.24)u"kJ/mol",
    "r5" => measurement(21.97, 2.3)u"kJ/mol",
    "r6" => measurement(-5.61, 0.55)u"kJ/mol",
    "biomass" => measurement(1080.17, 0.67)u"kJ/mol",
    "r7" => measurement(-2277.33, 1.49)u"kJ/mol",
    "r8" => measurement(1402.53, 1.23)u"kJ/mol",
    "r9" => measurement(0, 0)u"kJ/mol",
)

substrates = ["m1"]
targets = ["m6"]
cofactors = ["co1", "co2", "h"] #allow cofactors

os = build_OptStoic_model(m, substrates, targets, cofactors, Gurobi, energy_dc = formation_dG)
os = optimize_OptStoic(os)
if is_solved(os)
    os_solution = solution_OptStoic(os)
end
os_solution