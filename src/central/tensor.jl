"""
    CentralTensorMap(data::DenseVector, space::CentralTensorMapSpace)
    CentralTensorMap{T}(undef, codomain, domain)

Store the independent coefficients in a flat vector, without copying supplied
data. `block(U, a)` is a writable matrix in bulk sector `a`. `U[f, v, g]`
views its coefficients for a boundary splitting tree, junction, and bulk fusion
tree, with one degeneracy axis per external leg.
"""
struct CentralTensorMap{T,A<:DenseVector{T},S<:CentralTensorMapSpace}
    data::A
    space::S

    function CentralTensorMap(data::A, space::S) where {T,A<:DenseVector{T},S<:CentralTensorMapSpace}
        length(data) == dim(space) || throw(DimensionMismatch("expected $(dim(space)) coefficients, got $(length(data))"))
        return new{T,A,S}(data, space)
    end
end

CentralTensorMap(data::DenseVector, out::TensorSpace, input::TensorSpace) =
    CentralTensorMap(data, CentralTensorMapSpace(out, input))
Base.zeros(::Type{T}, S::CentralTensorMapSpace) where {T} = CentralTensorMap(zeros(T, dim(S)), S)

CentralTensorMap{T}(::UndefInitializer, S::CentralTensorMapSpace) where {T} =
    CentralTensorMap(Vector{T}(undef, dim(S)), S)
CentralTensorMap{T}(::UndefInitializer, out::TensorSpace, input::TensorSpace) where {T} =
    CentralTensorMap{T}(undef, CentralTensorMapSpace(out, input))

TensorKit.space(U::CentralTensorMap) = U.space
TensorKit.codomain(U::CentralTensorMap) = codomain(U.space)
TensorKit.domain(U::CentralTensorMap) = domain(U.space)
TensorKit.numout(U::CentralTensorMap) = TensorKit.numout(U.space)
TensorKit.numin(U::CentralTensorMap) = TensorKit.numin(U.space)
Base.eltype(::Type{<:CentralTensorMap{T}}) where {T} = T
TensorKit.scalartype(::Type{<:CentralTensorMap{T}}) where {T} = T
TensorKit.storagetype(::Type{<:CentralTensorMap{T,A}}) where {T,A} = A
TensorKit.blocksectors(U::CentralTensorMap) = blocksectors(U.space)
TensorKit.fusiontrees(U::CentralTensorMap) = fusiontrees(U.space)

_blockview(U::CentralTensorMap, B::CentralBlockLayout) = reshape(view(U.data, B.range), B.shape)

function TensorKit.block(U::CentralTensorMap, a::Sector)
    typeof(a) === sectortype(domain(U)) || throw(SectorMismatch("invalid bulk sector type"))
    layouts = U.space.layout
    haskey(layouts, a) && return _blockview(U, layouts[a])
    W = forget(a)
    m = sum(r -> blockdim(codomain(U), r) * dim(W, r), sectors(W); init=0)
    return reshape(view(U.data, 1:0), m, blockdim(domain(U), a))
end

function TensorKit.subblock(U::CentralTensorMap, (f, v, g)::Tuple{FusionTree,BoundaryBulkJunction,FusionTree})
    layouts = U.space.layout
    typeof(v.bulk) === sectortype(domain(U)) && haskey(layouts, v.bulk) ||
        throw(SectorMismatch("junction bulk sector is absent from this map"))
    B = layouts[v.bulk]
    haskey(B.rows, (f, v)) && haskey(B.columns, g) ||
        throw(SectorMismatch("trees and junction do not form a subblock of this map"))
    tile = view(_blockview(U, B), B.rows[(f, v)], B.columns[g])
    return sreshape(tile, (dims(codomain(U), f.uncoupled)..., dims(domain(U), g.uncoupled)...))
end

Base.getindex(U::CentralTensorMap, f::FusionTree, v::BoundaryBulkJunction, g::FusionTree) =
    TensorKit.subblock(U, (f, v, g))

function Base.setindex!(U::CentralTensorMap, A::AbstractArray, f::FusionTree, v::BoundaryBulkJunction, g::FusionTree)
    copy!(U[f, v, g], A)
    return U
end
