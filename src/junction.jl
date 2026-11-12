"""
    BoundaryBulkJunction(a, r, μ)

Label the standard projection from `forget(a)` onto its `μ`-th copy of the
simple boundary sector `r`. The copy index is one-based.
"""
struct BoundaryBulkJunction{I<:Sector,C<:DrinfeldCenter{I}}
    bulk::C
    boundary::I
    copyindex::Int

    @inline function BoundaryBulkJunction(a::C, r::I, μ::Integer) where {I<:Sector,C<:DrinfeldCenter{I}}
        @boundscheck begin
            copies = Base.OneTo(dim(forget(a), r))
            checkbounds(copies, μ)
        end
        return new{I,C}(a, r, μ)
    end
end

function Base.:(==)(v::BoundaryBulkJunction, w::BoundaryBulkJunction)
    return typeof(v) === typeof(w) && v.bulk == w.bulk &&
        v.boundary == w.boundary && v.copyindex == w.copyindex
end

function Base.hash(v::BoundaryBulkJunction, h::UInt)
    return hash((typeof(v), v.bulk, v.boundary, v.copyindex), h)
end

"""
    junctions(a)
    junctions(a, r)

Iterate the junction labels of a bulk sector in TensorKit's boundary-sector
order, then in increasing copy order. Pass `r` to select one boundary sector.
"""
function junctions(a::DrinfeldCenter)
    W = forget(a)
    return Iterators.flatten(_junctions(a, r, dim(W, r)) for r in sectors(W))
end

function junctions(a::DrinfeldCenter{I}, r::I) where {I<:Sector}
    return _junctions(a, r, dim(forget(a), r))
end

function _junctions(a::DrinfeldCenter{I}, r::I, n::Int) where {I<:Sector}
    return ((@inbounds BoundaryBulkJunction(a, r, μ)) for μ in Base.OneTo(n))
end

"""Materialize the junction as a one-output, one-input boundary tensor map."""
function TensorKit.TensorMap(v::BoundaryBulkJunction{I}) where {I}
    P = zeros(Vect[I](v.boundary => 1), forget(v.bulk))
    block(P, v.boundary)[1, v.copyindex] = 1
    return P
end
