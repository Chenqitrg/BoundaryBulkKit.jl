using Test
using TensorKit
using BoundaryBulkKit

include(joinpath(@__DIR__, "..", "examples", "support", "known_centers.jl"))
using .ExampleCenters: Z2Center, ElectricA4

@testset "Known center data in TensorKit bases" begin
    e, m = Z2Center(true, false), Z2Center(false, true)
    em = only(e ⊗ m)
    trivial, odd = Z2Irrep(0), Z2Irrep(1)

    @test base_sector(Z2Center) === base_sector(em) === Z2Irrep
    @test forget(em) == Vect[Z2Irrep](odd => 1)
    W = Vect[Z2Center](e => 2) ⊗ Vect[Z2Center](m => 3)
    f = only(fusiontrees(W, em))
    @test f.uncoupled == (e, m)
    @test f.coupled == em
    @test TensorKit.blockdim(W, em) == 6

    Ω = halfbraiding_symbol(em, odd, odd, odd, trivial, 1, 1)
    @test size(Ω) == (dim(forget(em), odd), dim(forget(em), odd))
    @test only(Ω) == -1
    @test only(halfbraiding_symbol(em, trivial, odd, odd, odd, 1, 1)) == 1
    vertex = forget_vertex(e, m, em, odd, trivial, odd, 1, 1)
    @test size(vertex) == (dim(forget(e), odd), dim(forget(m), trivial), dim(forget(em), odd))
    @test only(vertex) == 1
end

include("junction.jl")

include("centraltensormap.jl")
include("central_fusion.jl")

include("halfbraiding.jl")

include("center_adapters.jl")
