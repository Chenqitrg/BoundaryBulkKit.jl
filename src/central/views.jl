# Like TensorKit's blocks and subblocks, iterate key => view pairs and also
# allow direct lookup by key. Neither operation copies the coefficients.
struct CentralViews{F,U<:CentralTensorMap,K}
    lookup::F
    tensor::U
    labels::K
end

TensorKit.blocks(U::CentralTensorMap) = CentralViews(block, U, blocksectors(U))
TensorKit.subblocks(U::CentralTensorMap) = CentralViews(TensorKit.subblock, U, fusiontrees(U))

Base.getindex(V::CentralViews, key) = V.lookup(V.tensor, key)
Base.keys(V::CentralViews) = V.labels
Base.IteratorSize(::Type{<:CentralViews{F,U,K}}) where {F,U,K} = Base.IteratorSize(K)
Base.IteratorEltype(::Type{<:CentralViews}) = Base.EltypeUnknown()
Base.length(V::CentralViews) = length(V.labels)

function Base.iterate(V::CentralViews, state...)
    next = iterate(V.labels, state...)
    isnothing(next) && return nothing
    key, newstate = next
    return (key => V[key]), newstate
end
