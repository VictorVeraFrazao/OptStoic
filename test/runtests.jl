using OptStoic
using Test

@testset "OptStoic.jl test suite" begin
    include("thermodynamics.jl")
    include("analysis.jl")
    include("build_model.jl")
end
