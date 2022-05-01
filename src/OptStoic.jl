module OptStoic

using eQuilibrator
using Unitful, Measurements
using ProgressMeter
using Suppressor
using JuMP
using COBREXA
using OrderedCollections
using DocStringExtensions
using LinearAlgebra
using SparseArrays

include("thermodynamics.jl")
include("preprocessing.jl")
include("modify_model.jl")
include("build_model.jl")
include("Analysis.jl")

# export everything that isn't prefixed with _ (inspired by JuMP.jl and COBREXA.jl)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end
