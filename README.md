# OptStoic
This package provides the tool for the LP/MILP procedure OptStoic and MinFlux in Julia.
## Description
OptStoic is a two-step MP/MILP procedure for *de novo* metabolic pathway construction, the first step being an maximization of the stoichiometric coefficients of target molecules in the given overall biochemical reaction (OptStoic) and the second being a minimization of flux of the reactions in order to create a shortest corresponding pathway (MinFlux). Both steps are, inter alia, constrained by energy and mass balance.
## Dependencies
This package heavily depends on the [COBREXA package](https://github.com/LCSB-BioCore/COBREXA.jl) to handle models in Julia as well as [eQuilibrator package](https://github.com/stelmo/Equilibrator.jl) for Component Contribution methods required by the energy constraints. Currently, the package only is compatible with well-curated SBML models, provided by [BiGG database](bigg.ucsd.edu/) and the Gurobi optimizer. It is planned to extend the usage onto large-scale databases and include support for other LP solvers.
## Example usage
The package can be accessed by cloning the repository. The package can be activated as follows (dockerfiles will be added eventually):
```
using Pkg
Pkg.activate("path/to/repository/OptStoic/")
```
Then load the following packages, after installing them via the packaging environment (**] add PACKAGE**):
```
using OptStoic # Loads package
using COBREXA # Loads required metabolic model structures
using Gurobi # Loads solver. Keep in mind, that Gurobi is not open-source. Other feasible solvers may be tested
using eQuilibrator # Loads component contribution tools
using Unitful, Measurements # Loads physical unit support
```
We choose the basic E.coli core model and download it from the BiGG database:
```
download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")
model = load_model(StandardModel, "e_coli_core.xml")
```
The pathway assembly requires thermodynamic information that can be collected by using component contribution. The eQuilibrator package 
## Current developmental status
Work in progress, the code is partially stable. Basic workflow with BiGG models is possible, e.g. simple OptStoic+MinFlux pathway assembly (see Example usage).
## Reference
Please cite [Maranas, Costas. (2015). Designing overall stoichiometric conversions and intervening metabolic reactions. Scientific Reports. 5. 16009. 10.1038/srep16009. ](https://www.researchgate.net/publication/283979269_Designing_overall_stoichiometric_conversions_and_intervening_metabolic_reactions) and [Ng, C.Y., Wang, L., Chowdhury, A. et al. Pareto Optimality Explanation of the Glycolytic Alternatives in Nature. Sci Rep 9, 2633 (2019). https://doi.org/10.1038/s41598-019-38836-9](https://www.nature.com/articles/s41598-019-38836-9).
