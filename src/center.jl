"""
    DrinfeldCenter{I<:Sector} <: Sector

Supertype for the simple labels of the full Drinfeld center of boundary sector
type `I`. For a concrete center type `C <: DrinfeldCenter{I}`, `values(C)`
enumerates its simple objects, and `TensorKitSectors.ismodular(C)` returns
`true` by the center contract. This trait does not validate supplied symbols.

Implement TensorKit's sector interface for `C` (including bulk fusion and
braiding), then implement [`forget_Nsymbol`](@ref), [`halfbraiding_symbol`](@ref) and
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
    forget_Nsymbol(x::DrinfeldCenter{I}, a::I) where {I<:Sector}

Return the nonnegative integer multiplicity of boundary sector `a` in the
forgetful image of `x`, including zero for an absent sector. Copies are
numbered from one in the junction basis shared by both symbol methods.
"""
function forget_Nsymbol end

"""
    halfbraiding_symbol(a::C, i::I, r::I, s::I, j::I, μ::Int, ν::Int) -> AbstractMatrix{<:Number}

Implement for `I <: Sector` and `C <: DrinfeldCenter{I}`. Write `m(a,r) = forget_Nsymbol(a,r)`.
For the crossing `forget(a) ⊗ i → i ⊗ forget(a)` at total boundary charge `j`, the output
tree splits `j → i ⊗ r` in channel `1 ≤ μ ≤ Nsymbol(i,r,j)`. The input
tree splits `j → s ⊗ i` in channel `1 ≤ ν ≤ Nsymbol(s,i,j)`.

Return a matrix of size `(m(a,r), m(a,s))`. Its entry `[p,q]` maps
input junction copy `q` of `s` to output copy `p` of `r`. Even a single
phase occupies a 1×1 matrix. The caller supplies admissible tree channels.
"""
function halfbraiding_symbol end

"""
    forget_vertex(a::C, b::C, c::C, r::I, s::I, t::I, λ::Int, μ::Int) -> AbstractArray{<:Number,3}

Implement for `I <: Sector` and `C <: DrinfeldCenter{I}`. Expand the bulk
splitting vertex `c → a ⊗ b` (channel `1 ≤ λ ≤ Nsymbol(a,b,c)`) at the
boundary vertex `t → r ⊗ s` (channel `1 ≤ μ ≤ Nsymbol(r,s,t)`).

With `m = forget_Nsymbol`, return shape `(m(a,r), m(b,s), m(c,t))`.
Entry `[p,q,k]` goes from input copy `k` of `t` to output copies `p,q`
of `r,s`. The caller supplies admissible tree channels.
"""
function forget_vertex end

"""
    forgetvertexscalartype(::Type{C}) where {C<:DrinfeldCenter}

Scalar type for coefficients of the forgetful splitting vertices. As with
`TensorKit.fusionscalartype`, the default uses the unit vertex. Specialize
this method when other vertices need a wider type or the unit vertex is
expensive to construct.
"""
function forgetvertexscalartype(::Type{C})::Type{<:Number} where {I<:Sector,C<:DrinfeldCenter{I}}
    x, a = unit(C), unit(I)
    return eltype(forget_vertex(x, x, x, a, a, a, 1, 1))
end

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
