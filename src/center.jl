"""
    DrinfeldCenter{I<:Sector} <: Sector

Supertype for the simple labels of the full Drinfeld center of boundary sector
type `I`. For a concrete center type `C <: DrinfeldCenter{I}`, `values(C)`
enumerates its simple objects, and `TensorKitSectors.ismodular(C)` returns
`true` by the center contract. This trait does not validate supplied symbols.

Implement TensorKit's sector interface for `C` (including bulk fusion and
braiding), then implement [`forget`](@ref), [`halfbraiding_symbol`](@ref) and
[`forget_vertex`](@ref) in the same junction and fusion-tree bases.
"""
abstract type DrinfeldCenter{I<:Sector} <: Sector end

TKS.ismodular(::Type{<:DrinfeldCenter}; kwargs...) = true

"""
    base_sector(::Type{C}) -> Type{I}
    base_sector(a::C) -> Type{I}

For `C <: DrinfeldCenter{I}`, return the boundary sector type `I`.
"""
function base_sector(::Type{<:DrinfeldCenter{I}})::Type{I} where {I<:Sector}
    return I
end
function base_sector(::DrinfeldCenter{I})::Type{I} where {I<:Sector}
    return I
end

"""
    forget(a::C) -> ElementarySpace

Implement for `C <: DrinfeldCenter{I}`. The result `X` is one boundary leg
with `sectortype(X) === I`. The built-in centers use `Vect[I]`, which is
`ComplexSpace` for `Trivial` and a concrete `GradedSpace{I}` otherwise.
`dim(X, r)` is the multiplicity of `r::I`; the integers `1:dim(X, r)`
index its junction basis. Both symbol methods use this copy order.
"""
function forget end

"""
    halfbraiding_symbol(a::C, i::I, r::I, s::I, j::I, μ::Int, ν::Int) -> AbstractMatrix{<:Number}

Implement for `I <: Sector` and `C <: DrinfeldCenter{I}`. Write `X = forget(a)`.
For the crossing `X ⊗ i → i ⊗ X` at total boundary charge `j`, the output
tree splits `j → i ⊗ r` in channel `1 ≤ μ ≤ Nsymbol(i,r,j)`. The input
tree splits `j → s ⊗ i` in channel `1 ≤ ν ≤ Nsymbol(s,i,j)`.

Return a matrix of size `(dim(X,r), dim(X,s))`. Its entry `[p,q]` maps
input junction copy `q` of `s` to output copy `p` of `r`. Even a single
phase occupies a 1×1 matrix. The caller supplies admissible tree channels.
"""
function halfbraiding_symbol end

"""
    forget_vertex(a::C, b::C, c::C, r::I, s::I, t::I, λ::Int, μ::Int) -> AbstractArray{<:Number,3}

Implement for `I <: Sector` and `C <: DrinfeldCenter{I}`. Expand the bulk
splitting vertex `c → a ⊗ b` (channel `1 ≤ λ ≤ Nsymbol(a,b,c)`) at the
boundary vertex `t → r ⊗ s` (channel `1 ≤ μ ≤ Nsymbol(r,s,t)`).

Return shape `(dim(forget(a),r), dim(forget(b),s), dim(forget(c),t))`.
Entry `[p,q,k]` goes from input copy `k` of `t` to output copies `p,q`
of `r,s`. The caller supplies admissible tree channels.
"""
function forget_vertex end

# The boundary irrep type already carries its group parameter.
abstract type GroupDoubleSector{I<:TKS.AbstractIrrep} <: DrinfeldCenter{I} end

"""
    QuantumDouble{G}

Family of group quantum doubles for the TensorKit group type `G`.
The concrete boundary irrep type can be recovered with `base_sector`.
For example, `Z3QuantumDouble <: QuantumDouble{ℤ{3}}` and its boundary
sector type is `Z3Irrep`.
"""
const QuantumDouble{G<:TKS.Group} = GroupDoubleSector{I} where {I<:TKS.AbstractIrrep{G}}

"""
    take_center(::Type{I}) -> Type{<:DrinfeldCenter{I}}

Return the registered center type for the boundary sector type `I <: Sector`.
Built-in choices are the trivial center for `Trivial`, `ZNQuantumDouble{N}`
for `ZNIrrep{N}`, and `ModularDouble{I}` for modular boundaries. Other
boundaries require a registration; the generic method throws `ArgumentError`
if the modularity check returns false.
"""
function take_center end
