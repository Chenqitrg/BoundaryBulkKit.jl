# Within each boundary sector, order by bulk sector, bulk copy, then junction.
function _forget_layout(A::GradedSpace{C}) where {I<:Sector,C<:DrinfeldCenter{I}}
    totals = Dict{I,Int}()
    summands = map(collect(sectors(A))) do a
        X = forget(a)
        copies = dim(A, a)
        firstcopy = Dict{I,UnitRange{Int}}()
        for r in sectors(X)
            offset = get(totals, r, 0)
            n = dim(X, r)
            firstcopy[r] = (offset + 1):(offset + n)
            totals[r] = offset + copies * n
        end
        return (; a, copies, firstcopy)
    end
    W = Vect[I](totals)
    return (isdual(A) ? flip(W) : W), summands
end

"""
    forget(A::GradedSpace{C}) -> ElementarySpace, with C <: DrinfeldCenter{I}

Underlying boundary space of a center leg. Within each boundary sector, copies
are ordered by `sectors(A)`, then by copy in `A`, then by junction in the
simple object's fixed basis. The result has sector type `I` and the same dual
orientation as `A`.
"""
function forget(A::GradedSpace{C})::Vect[I] where {I<:Sector,C<:DrinfeldCenter{I}}
    return first(_forget_layout(A))
end

"""
    HalfBraiding(A::GradedSpace{C}, V::ElementarySpace; T=ComplexF64) -> TensorMap
    HalfBraiding(a::C, V::ElementarySpace; T=ComplexF64) -> TensorMap

Materialize the boundary TensorMap `V ⊗ forget(A) ← forget(A) ⊗ V`, for a
center-graded leg `A` and a boundary leg `V`. Each simple half-braiding
symbol acts on junction copies and preserves both external multiplicity
indices. For `C <: DrinfeldCenter{I}`, the boundary leg must have sector type
`I`. A simple label `a` represents one copy. The result has two output and two
input legs. Choose a numeric `T` that can represent the crossing coefficients.
"""
function HalfBraiding(A::GradedSpace{C}, V::ElementarySpace; T::Type{<:Number}=ComplexF64)::TensorMap where {I<:Sector,C<:DrinfeldCenter{I}}
    sectortype(V) === I || throw(SectorMismatch("boundary leg must have sector type $I"))
    X, summands = _forget_layout(A)
    return _assemble_halfbraiding(T, X, V, summands)
end

function HalfBraiding(a::DrinfeldCenter{I}, V::ElementarySpace; T::Type{<:Number}=ComplexF64)::TensorMap where {I<:Sector}
    sectortype(V) === I || throw(SectorMismatch("boundary leg must have sector type $I"))
    X = forget(a)
    firstcopy = Dict(r => 1:dim(X, r) for r in sectors(X))
    return _assemble_halfbraiding(T, X, V, ((; a, copies=1, firstcopy),))
end

function _assemble_halfbraiding(::Type{T}, X::ElementarySpace, V::ElementarySpace, summands) where {T}
    Γ = zeros(T, V ⊗ X ← X ⊗ V)
    for (f, g) in fusiontrees(Γ)
        i, r = f.uncoupled
        s, k = g.uncoupled
        i == k || continue
        B = Γ[f, g]
        for part in summands
            rows = get(part.firstcopy, r, 1:0)
            cols = get(part.firstcopy, s, 1:0)
            (isempty(rows) || isempty(cols)) && continue
            Ω = halfbraiding_symbol(part.a, i, r, s, f.coupled, only(f.vertices), only(g.vertices))
            size(Ω) == (length(rows), length(cols)) ||
                throw(DimensionMismatch("half-braiding symbol has size $(size(Ω)); expected $((length(rows), length(cols)))"))
            for p in 0:(part.copies - 1)
                out = rows .+ p * length(rows)
                input = cols .+ p * length(cols)
                for b in axes(B, 1)
                    B[b, out, input, b] .= Ω
                end
            end
        end
    end
    return Γ
end
