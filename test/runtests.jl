using Test
using TensorKit
using BoundaryBulkKit

# Known data for the untwisted Z₂ double in the charge-flux basis.
struct Z2Center <: DrinfeldCenter{Z2Irrep}
    charge::Bool
    flux::Bool
end

const Z2_LABELS = tuple((Z2Center(q, m) for q in (false, true) for m in (false, true))...)
const Z2Labels = typeof(values(Z2Center))
Base.length(::Z2Labels) = length(Z2_LABELS)
Base.getindex(::Z2Labels, i::Int) = Z2_LABELS[i]
Base.iterate(::Z2Labels, i::Int=1) = iterate(Z2_LABELS, i)

TensorKit.unit(::Type{Z2Center}) = Z2Center(false, false)
TensorKit.dual(a::Z2Center) = a
TensorKit.FusionStyle(::Type{Z2Center}) = UniqueFusion()
TensorKit.:⊗(a::Z2Center, b::Z2Center) = (Z2Center(a.charge ⊻ b.charge, a.flux ⊻ b.flux),)
TensorKit.Nsymbol(a::Z2Center, b::Z2Center, c::Z2Center) = c == only(a ⊗ b)
Base.isless(a::Z2Center, b::Z2Center) = isless((a.charge, a.flux), (b.charge, b.flux))

BoundaryBulkKit.forget(a::Z2Center) = Vect[Z2Irrep](Z2Irrep(a.charge) => 1)

function BoundaryBulkKit.halfbraiding_symbol(a::Z2Center, i::Z2Irrep, r::Z2Irrep,
        s::Z2Irrep, j::Z2Irrep, μ::Int, ν::Int)
    @assert r == s == Z2Irrep(a.charge)
    @assert μ == ν == 1 && j == only(i ⊗ r)
    return fill(a.flux && i == Z2Irrep(1) ? -1 : 1, 1, 1)
end

function BoundaryBulkKit.forget_vertex(a::Z2Center, b::Z2Center, c::Z2Center,
        r::Z2Irrep, s::Z2Irrep, t::Z2Irrep, λ::Int, μ::Int)
    @assert λ == μ == 1 && c == only(a ⊗ b)
    @assert (r, s, t) == (Z2Irrep(a.charge), Z2Irrep(b.charge), Z2Irrep(c.charge))
    return ones(Int, 1, 1, 1)
end

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
