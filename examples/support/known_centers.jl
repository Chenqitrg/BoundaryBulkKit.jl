module ExampleCenters

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

# Pure electric objects in Z(Rep(A₄)) inherit their fusion from Rep(A₄).
struct ElectricA4 <: DrinfeldCenter{A4Irrep}
    rep::A4Irrep
end
Base.length(::typeof(values(ElectricA4))) = length(values(A4Irrep))
Base.getindex(::typeof(values(ElectricA4)), i::Int) = ElectricA4(values(A4Irrep)[i])
Base.iterate(v::typeof(values(ElectricA4)), i::Int=1) = i <= length(v) ? (v[i], i + 1) : nothing
Base.isless(a::ElectricA4, b::ElectricA4) = isless(a.rep, b.rep)
TensorKit.unit(::Type{ElectricA4}) = ElectricA4(unit(A4Irrep))
TensorKit.dual(a::ElectricA4) = ElectricA4(dual(a.rep))
TensorKit.FusionStyle(::Type{ElectricA4}) = FusionStyle(A4Irrep)
TensorKit.:⊗(a::ElectricA4, b::ElectricA4) = (ElectricA4(r) for r in a.rep ⊗ b.rep)
TensorKit.Nsymbol(a::ElectricA4, b::ElectricA4, c::ElectricA4) = Nsymbol(a.rep, b.rep, c.rep)
BoundaryBulkKit.forget(a::ElectricA4) = Vect[A4Irrep](a.rep => 1)

end
