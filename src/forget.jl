"""Construct the boundary space of a simple center object from its junction multiplicities."""
function forget(x::DrinfeldCenter{I})::Vect[I] where {I<:Sector}
    return Vect[I](a => forget_Nsymbol(x, a) for a in _boundary_labels(x))
end

function _boundary_labels(::DrinfeldCenter{I})::TKS.SectorSet{I} where {I<:Sector}
    return TKS.SectorSet{I}(values(I))
end

const _JunctionRanges{C,I} = Dict{C,Dict{I,UnitRange{Int}}}

# Within each boundary sector, order by bulk sector, bulk copy, then junction.
function _forget_layout(A::GradedSpace{C})::Tuple{Vect[I],Vector{<:NamedTuple}} where {I<:Sector,C<:DrinfeldCenter{I}}
    totals = Dict{I,Int}()
    summands = map(collect(sectors(A))) do a::C
        copies = dim(A, a)
        firstcopy = Dict{I,UnitRange{Int}}()
        for r in _boundary_labels(a)
            n = forget_Nsymbol(a, r)
            iszero(n) && continue
            n >= 0 || throw(ArgumentError("negative forget_Nsymbol for $a and $r"))
            offset = get(totals, r, 0)
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

# Matrix rows index all leaf junctions (first leg varies fastest); columns
# index junctions on the root. Contract adjacent vertices on their shared edge.
function _tree_coefficients(::Type{T}, f::FusionTree{C,N}, h::FusionTree{I,N},
                            ranges::NTuple{N,UnitRange{Int}})::Union{Nothing,Matrix{T}} where
                           {T<:Number,I<:Sector,C<:DrinfeldCenter{I},N}
    isempty(ranges) && return ones(T, 1, 1)
    left, r = first(f.uncoupled), first(h.uncoupled)
    n = length(first(ranges))
    coeff::Matrix{T} = [T(i == j) for i in 1:n, j in 1:n]
    bulk_roots = (f.innerlines..., f.coupled)
    boundary_roots = (h.innerlines..., h.coupled)
    for i in eachindex(f.vertices)
        root, t = bulk_roots[i], boundary_roots[i]
        copies = forget_Nsymbol(root, t)
        iszero(copies) && return nothing
        vertex = forget_vertex(left, f.uncoupled[i+1], root, r, h.uncoupled[i+1],
                               t, f.vertices[i], h.vertices[i])
        size(vertex) == (size(coeff, 2), length(ranges[i+1]), copies) ||
            throw(DimensionMismatch("inconsistent junction dimensions in forget_vertex"))
        coeff = reshape(coeff * reshape(vertex, size(vertex, 1), :), :, copies)
        left, r = root, t
    end
    return coeff
end

struct _TreeTerm{F<:FusionTree,N,T<:Number}
    tree::F
    ranges::NTuple{N,UnitRange{Int}}
    coeff::Matrix{T}
end

function _forget_product(::Type{T}, P::ProductSpace{S,N}, mapspace::AbstractTensorMap{<:Number,S})::Tuple{ProductSpace,Dict} where
                        {T<:Number,I<:Sector,C<:DrinfeldCenter{I},S<:GradedSpace{C},N}
    parts = map(_forget_layout, P)
    W = ProductSpace{Vect[I]}(map(first, parts))
    layouts = map(parts) do part::Tuple{ElementarySpace,Vector{<:NamedTuple}}
        _JunctionRanges{C,I}(s.a => s.firstcopy for s in last(part))
    end

    F = TensorKit.fusiontreetype(C, N)
    Term = _TreeTerm{TensorKit.fusiontreetype(I, N),N,T}
    images = Dict{F,Vector{Term}}()
    for x in blocksectors(mapspace), f in fusiontrees(P, x)
        terms = Term[]
        for a in sectors(forget(x)), h in fusiontrees(W, a)
            ranges = map(layouts, f.uncoupled, h.uncoupled) do layout::_JunctionRanges{C,I}, y::C, b::I
                get(layout[y], b, 1:0)
            end
            any(isempty, ranges) && continue
            coeff = _tree_coefficients(T, f, h, ranges)
            isnothing(coeff) || push!(terms, Term(h, ranges, coeff))
        end
        images[f] = terms
    end
    return W, images
end

function _add_tree_image!(dest::AbstractArray{<:Number}, src::AbstractArray{<:Number},
                          rows::Tuple{Vararg{UnitRange{Int}}}, cols::Tuple{Vararg{UnitRange{Int}}},
                          weights::AbstractMatrix{<:Number})::Nothing
    ranges = (rows..., cols...)
    if isempty(ranges)
        dest[] += src[] * only(weights)
        return
    end
    sizes = map(length, ranges)
    junctions = reshape(weights, sizes)
    for index in CartesianIndices(src)
        amplitude = src[index]
        iszero(amplitude) && continue
        slab = map(ranges, sizes, Tuple(index)) do range::UnitRange{Int}, n::Int, k::Int
            range .+ (k - 1) * n
        end
        view(dest, slab...) .+= amplitude .* junctions
    end
end

"""
    forget(T::AbstractTensorMap{E,S}) where {E,C<:DrinfeldCenter,S<:GradedSpace{C}}

Apply the forgetful functor to a center-graded map. Each boundary leg uses the
same copy order as `forget(V)`. Output splitting trees are expanded using
`forget_vertex`; input trees use the adjoint expansion.
"""
function forget(T::AbstractTensorMap{E,S,M,N})::TensorMap{<:Number,Vect[I],M,N} where
               {E<:Number,I<:Sector,C<:DrinfeldCenter{I},S<:GradedSpace{C},M,N}
    # Zero- and one-leg trees have only identity junction coefficients.
    coefftype = max(M, N) < 2 ? Int : promote_type(Int, forgetvertexscalartype(C))
    out, rows = _forget_product(coefftype, codomain(T), T)
    input, cols = _forget_product(coefftype, domain(T), T)
    result = zeros(promote_type(E, coefftype), out ← input)
    for (f, g) in fusiontrees(T), row in rows[f], col in cols[g]
        row.tree.coupled == col.tree.coupled || continue
        weights = row.coeff * col.coeff'
        _add_tree_image!(result[row.tree, col.tree], T[f, g], row.ranges, col.ranges, weights)
    end
    return result
end
