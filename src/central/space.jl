# Within a bulk sector the columns have the TensorKit fusion-tree order of the
# domain. Rows have the same order as the coefficients of a boundary state
# with space codomain ← forget(a): boundary sector, copy, then splitting tree.
struct CentralBlockLayout{R,F}
    rows::Dictionary{R,UnitRange{Int}}
    columns::Dictionary{F,UnitRange{Int}}
    shape::Dims{2}
    range::UnitRange{Int}
end

function _axisranges(width, ::Type{K}, labels) where {K}
    ranges = Dictionary{K,UnitRange{Int}}()
    offset = 0
    for label in labels
        n = width(label)
        insert!(ranges, label, (offset + 1):(offset + n))
        offset += n
    end
    return ranges, offset
end

function _centralblocks(out::ProductSpace, input::ProductSpace,
                        ::Type{I}, ::Type{C}) where {I<:Sector,C<:DrinfeldCenter{I}}
    Fout = TensorKit.fusiontreetype(I, length(out))
    Fin = TensorKit.fusiontreetype(C, length(input))
    R = Tuple{Fout,BoundaryBulkJunction{I,C}}
    layouts = Dictionary{C,CentralBlockLayout{R,Fin}}()
    offset = 0
    for a in blocksectors(input)
        W = forget(a)
        rowlabels = ((f, v) for r in blocksectors(ProductSpace(W))
                     for v in _junctions(a, r, dim(W, r)) for f in fusiontrees(out, r))
        rows, m = _axisranges(R, rowlabels) do (f, v)
            dim(out, f.uncoupled)
        end
        iszero(m) && continue
        columns, n = _axisranges(Fin, fusiontrees(input, a)) do g
            dim(input, g.uncoupled)
        end
        range = (offset + 1):(offset + m * n)
        insert!(layouts, a, CentralBlockLayout(rows, columns, (m, n), range))
        offset += m * n
    end
    return layouts
end

"""
    CentralTensorMapSpace(codomain, domain)

Space of maps from bulk legs to boundary legs. Both arguments may be elementary
spaces or TensorKit product spaces. Construct it once to share the block layout
between maps with the same spaces. The leading type parameters are the boundary
sector `I` of the codomain and the bulk sector `C <: DrinfeldCenter{I}` of the domain.
`N` counts the bulk legs, including zero for the tensor unit.
"""
struct CentralTensorMapSpace{I<:Sector,C<:DrinfeldCenter{I},N,O<:ProductSpace,B<:ElementarySpace,
                             L<:Dictionary{C,<:CentralBlockLayout}}
    codomain::O
    domain::ProductSpace{B,N}
    layout::L

    function CentralTensorMapSpace(codomain::TensorSpace, domain::TensorSpace)
        out, input = ProductSpace(codomain), ProductSpace(domain)
        I, C = sectortype(out), sectortype(input)
        C <: DrinfeldCenter{I} || throw(SectorMismatch("bulk sector $C must have boundary sector $I"))
        layout = _centralblocks(out, input, I, C)
        return new{I,C,length(input),typeof(out),spacetype(input),typeof(layout)}(out, input, layout)
    end
end

TensorKit.codomain(S::CentralTensorMapSpace) = S.codomain
TensorKit.domain(S::CentralTensorMapSpace) = S.domain
TensorKit.numout(S::CentralTensorMapSpace) = length(codomain(S))
TensorKit.numin(S::CentralTensorMapSpace) = length(domain(S))
TensorKit.dim(S::CentralTensorMapSpace) = sum(B -> length(B.range), S.layout; init=0)
TensorKit.blocksectors(S::CentralTensorMapSpace) = keys(S.layout)
Base.:(==)(S::CentralTensorMapSpace, T::CentralTensorMapSpace) =
    codomain(S) == codomain(T) && domain(S) == domain(T)
Base.hash(S::CentralTensorMapSpace, h::UInt) = hash((codomain(S), domain(S)), h)

# Keep the Cartesian product of row and column labels lazy, even for large blocks.
struct CentralTreeTriples{F,V,G,L}
    labels::L
    count::Int
end
Base.eltype(::Type{<:CentralTreeTriples{F,V,G}}) where {F,V,G} = Tuple{F,V,G}
Base.length(K::CentralTreeTriples) = K.count
Base.iterate(K::CentralTreeTriples, state...) = iterate(K.labels, state...)

TensorKit.fusiontrees(S::CentralTensorMapSpace) = _treetriples(S.layout)
function _treetriples(layouts::Dictionary{C,CentralBlockLayout{Tuple{F,V},G}}) where {C,F,V,G}
    labels = ((f, v, g) for B in layouts for g in keys(B.columns) for (f, v) in keys(B.rows))
    n = sum(B -> length(B.rows) * length(B.columns), layouts; init=0)
    return CentralTreeTriples{F,V,G,typeof(labels)}(labels, n)
end
