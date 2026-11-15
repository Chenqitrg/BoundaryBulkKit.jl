import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
Pkg.instantiate(; allow_autoprecomp=false, io=devnull)

using TensorKit, BoundaryBulkKit
