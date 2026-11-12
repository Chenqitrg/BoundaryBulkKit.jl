# VectorInterface operations used by matrix-free Krylov methods.
TensorKit.zerovector(x::CentralTensorMap, ::Type{T}) where {T<:Number} = zeros(T, space(x))
TensorKit.scale(x::CentralTensorMap, α::Number) = CentralTensorMap(TensorKit.scale(x.data, α), space(x))
function TensorKit.scale!!(x::CentralTensorMap, α::Number)
    data = TensorKit.scale!!(x.data, α)
    return data === x.data ? x : CentralTensorMap(data, space(x))
end
function TensorKit.add!!(y::CentralTensorMap, x::CentralTensorMap, α::Number, β::Number)
    space(y) == space(x) || throw(TensorKit.SpaceMismatch("central vector spaces differ"))
    data = TensorKit.add!!(y.data, x.data, α, β)
    return data === y.data ? y : CentralTensorMap(data, space(y))
end
function TensorKit.scale!!(y::CentralTensorMap, x::CentralTensorMap, α::Number)
    space(y) == space(x) || throw(TensorKit.SpaceMismatch("central vector spaces differ"))
    data = TensorKit.scale!!(y.data, x.data, α)
    return data === y.data ? y : CentralTensorMap(data, space(y))
end

# The boundary categorical trace weights each row by its coupled boundary dimension.
function TensorKit.inner(x::CentralTensorMap, y::CentralTensorMap)
    space(x) == space(y) || throw(TensorKit.SpaceMismatch("central vector spaces differ"))
    result = zero(promote_type(eltype(x), eltype(y)))
    for a in blocksectors(x)
        X, Y = block(x,a), block(y,a)
        for ((_, v), rows) in pairs(space(x).layout[a].rows)
            result += dim(v.boundary) * TensorKit.dot(view(X,rows,:), view(Y,rows,:))
        end
    end
    return result
end
TensorKit.norm(x::CentralTensorMap) = sqrt(real(TensorKit.inner(x,x)))
Base.:*(α::Number, x::CentralTensorMap) = TensorKit.scale(x, α)
Base.isapprox(x::CentralTensorMap, y::CentralTensorMap; kwargs...) =
    space(x) == space(y) && isapprox(x.data,y.data; kwargs...)

function Base.:*(U::CentralTensorMap, D::TensorKit.DiagonalTensorMap)
    domain(U) == codomain(D) || throw(TensorKit.SpaceMismatch("incompatible central and diagonal factors"))
    result = CentralTensorMap{promote_type(eltype(U),eltype(D))}(undef, space(U))
    for a in blocksectors(U)
        block(result,a) .= block(U,a) .* transpose(block(D,a).diag)
    end
    return result
end
