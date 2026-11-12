"""
    DrinfeldCenter{I} <: Sector

Supertype of simple bulk sectors in the Drinfeld center of the boundary sector
type `I`. Concrete subtypes implement TensorKit's ordinary sector interface in
the bulk and supply the three bulk–boundary operations below in a fixed basis.
"""
abstract type DrinfeldCenter{I<:Sector} <: Sector end

"""
    base_sector(C)
    base_sector(a)

Boundary sector type associated with a concrete center sector type or label.
"""
base_sector(::Type{<:DrinfeldCenter{I}}) where {I} = I
base_sector(a::DrinfeldCenter) = base_sector(typeof(a))

"""
    forget(a::DrinfeldCenter{I})

Return the underlying single boundary leg as a TensorKit space graded by `I`.
The copies of a boundary sector `r`, indexed from 1 to `dim(forget(a), r)`,
label the fixed junction basis of projections from the underlying object to `r`.
Both symbols below use this same ordering.
"""
function forget end

"""
    halfbraiding_symbol(a, i, r, s, j, μ, ν)

Matrix of the crossing `forget(a) ⊗ i → i ⊗ forget(a)` at simple boundary
sector `i` and total boundary sector `j`. In TensorKit's splitting-tree basis
the output vertex is `j → i ⊗ r` (channel `μ`) and the input tree is
`j → s ⊗ i` (channel `ν`).

Return an `AbstractMatrix` of size `(dim(forget(a), r), dim(forget(a), s))`.
The row selects a copy of `r` on the output, the column a copy of `s` on the
input. Channel indices are one-based. Only admissible vertices are queried.
"""
function halfbraiding_symbol end

"""
    forget_vertex(a, b, c, r, s, t, λ, μ)

Boundary coefficients of the bulk splitting vertex `c → a ⊗ b` (channel `λ`)
on the boundary splitting vertex `t → r ⊗ s` (channel `μ`).

Return a three-axis array indexed by the copies of `r` in `forget(a)`, `s` in
`forget(b)`, and `t` in `forget(c)`, in that order. Channel indices are
one-based. Only admissible vertices are queried.
"""
function forget_vertex end
