import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
Pkg.instantiate(; allow_autoprecomp=false, io=devnull)

using TensorKit, BoundaryBulkKit
if !isdefined(@__MODULE__, :ExampleCenters)
    include(joinpath(@__DIR__, "known_centers.jl"))
end
using .ExampleCenters: Z2Center, ElectricA4
