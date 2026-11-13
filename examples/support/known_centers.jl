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
TensorKit.BraidingStyle(::Type{Z2Center}) = Anyonic()
TensorKit.Fsymbol(a::Z2Center,b::Z2Center,c::Z2Center,d::Z2Center,e::Z2Center,f::Z2Center) =
    Nsymbol(a,b,e) && Nsymbol(e,c,d) && Nsymbol(b,c,f) && Nsymbol(a,f,d)
TensorKit.Rsymbol(a::Z2Center,b::Z2Center,c::Z2Center) =
    Nsymbol(a,b,c) * (a.flux && b.charge ? -1 : 1)
TensorKit.findindex(::Z2Labels, a::Z2Center) = 1 + 2 * a.charge + a.flux
TensorKit.FusionStyle(::Type{Z2Center}) = UniqueFusion()
TensorKit.:⊗(a::Z2Center, b::Z2Center) = (Z2Center(a.charge ⊻ b.charge, a.flux ⊻ b.flux),)
TensorKit.Nsymbol(a::Z2Center, b::Z2Center, c::Z2Center) = c == only(a ⊗ b)
Base.isless(a::Z2Center, b::Z2Center) = isless((a.charge, a.flux), (b.charge, b.flux))

BoundaryBulkKit.forget(a::Z2Center)::GradedSpace{Z2Irrep} = Vect[Z2Irrep](Z2Irrep(a.charge) => 1)

function BoundaryBulkKit.halfbraiding_symbol(a::Z2Center, i::Z2Irrep, r::Z2Irrep,
        s::Z2Irrep, j::Z2Irrep, μ::Int, ν::Int)::Matrix{Int}
    @assert r == s == Z2Irrep(a.charge)
    @assert μ == ν == 1 && j == only(i ⊗ r)
    return fill(a.flux && i == Z2Irrep(1) ? -1 : 1, 1, 1)
end

function BoundaryBulkKit.forget_vertex(a::Z2Center, b::Z2Center, c::Z2Center,
        r::Z2Irrep, s::Z2Irrep, t::Z2Irrep, λ::Int, μ::Int)::Array{Int,3}
    @assert λ == μ == 1 && c == only(a ⊗ b)
    @assert (r, s, t) == (Z2Irrep(a.charge), Z2Irrep(b.charge), Z2Irrep(c.charge))
    return ones(Int, 1, 1, 1)
end

end
