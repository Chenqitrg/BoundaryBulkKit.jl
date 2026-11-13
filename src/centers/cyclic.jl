"""Charge and flux of the untwisted quantum double of ℤₙ."""
struct ZNQuantumDouble{N} <: GroupDoubleSector{ZNIrrep{N}}
    charge::Int
    flux::Int
    function ZNQuantumDouble{N}(q::Integer, m::Integer) where {N}
        return new{N}(mod(q, N), mod(m, N))
    end
end

const Z2QuantumDouble = ZNQuantumDouble{2}
const Z3QuantumDouble = ZNQuantumDouble{3}
const 𝒟ZN = ZNQuantumDouble
const 𝒟Z2 = Z2QuantumDouble
const 𝒟Z3 = Z3QuantumDouble

function take_center(::Type{ZNIrrep{N}})::Type{ZNQuantumDouble{N}} where {N}
    return ZNQuantumDouble{N}
end
TensorKit.unit(::Type{C}) where {C<:ZNQuantumDouble} = C(0, 0)
TensorKit.dual(a::C) where {C<:ZNQuantumDouble} = C(-a.charge, -a.flux)
TensorKit.:⊗(a::C, b::C) where {C<:ZNQuantumDouble} = (C(a.charge+b.charge, a.flux+b.flux),)
TensorKit.Nsymbol(a::C,b::C,c::C) where {C<:ZNQuantumDouble} = c == only(a ⊗ b)
TensorKit.FusionStyle(::Type{<:ZNQuantumDouble}) = TensorKit.UniqueFusion()
TensorKit.BraidingStyle(::Type{<:ZNQuantumDouble}) = TensorKit.Anyonic()
function TensorKit.Fsymbol(a::C,b::C,c::C,d::C,e::C,f::C) where {C<:ZNQuantumDouble}
    return e == only(a ⊗ b) && d == only(e ⊗ c) && f == only(b ⊗ c) && d == only(a ⊗ f)
end
TensorKit.Rsymbol(a::C,b::C,c::C) where {N,C<:ZNQuantumDouble{N}} =
    TensorKit.Nsymbol(a,b,c) * cispi(2 * a.flux * b.charge / N)
Base.isless(a::C,b::C) where {C<:ZNQuantumDouble} = isless((a.charge,a.flux),(b.charge,b.flux))
Base.length(::TKS.SectorValues{ZNQuantumDouble{N}}) where {N} = N^2
function Base.getindex(v::TKS.SectorValues{ZNQuantumDouble{N}}, i::Int) where {N}
    1 <= i <= length(v) || throw(BoundsError(v,i))
    return ZNQuantumDouble{N}(divrem(i-1,N)...)
end
Base.iterate(v::TKS.SectorValues{C}, i::Int=1) where {C<:ZNQuantumDouble} =
    i <= length(v) ? (v[i], i+1) : nothing
TKS.findindex(::TKS.SectorValues{ZNQuantumDouble{N}}, a::ZNQuantumDouble{N}) where {N} = 1+a.charge*N+a.flux

function forget(a::ZNQuantumDouble{N})::GradedSpace{ZNIrrep{N}} where {N}
    return Vect[ZNIrrep{N}](ZNIrrep{N}(a.charge)=>1)
end
function halfbraiding_symbol(a::ZNQuantumDouble{N}, i::I,r::I,s::I,j::I,μ::Int,ν::Int)::Matrix{ComplexF64} where {N,I<:ZNIrrep{N}}
    return fill(cispi(2 * a.flux * i.n / N), 1, 1)
end
function forget_vertex(a::C,b::C,c::C,r::I,s::I,t::I,λ::Int,μ::Int)::Array{Int,3} where {N,C<:ZNQuantumDouble{N},I<:ZNIrrep{N}}
    return ones(Int,1,1,1)
end
