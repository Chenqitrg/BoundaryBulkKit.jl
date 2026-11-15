"""
    ModularDouble(left::I, right::I) where {I<:Sector}

Construct a simple object of the double from two sectors of the same modular
category. The right factor carries the reverse braiding.
"""
struct ModularDouble{I<:Sector} <: DrinfeldCenter{I}
    label::ProductSector{Tuple{I,TimeReversed{I}}}

    function ModularDouble{I}(left::I, right::I)::ModularDouble{I} where {I<:Sector}
        TKS.ismodular(I) || throw(ArgumentError("ModularDouble requires a modular boundary sector, got $I"))
        return new{I}(ProductSector((left, TimeReversed{I}(right))))
    end
end
function ModularDouble(left::I, right::I)::ModularDouble{I} where {I<:Sector}
    return ModularDouble{I}(left, right)
end

function _from_product(label::ProductSector{Tuple{I,TimeReversed{I}}})::ModularDouble{I} where {I<:Sector}
    left, right = label.sectors
    return ModularDouble(left, right.a)
end

function Base.show(io::IO, x::ModularDouble)::Nothing
    left, right = x.label.sectors
    print(io, "ModularDouble(")
    show(io, left)
    print(io, ", ")
    show(io, right.a)
    print(io, ")")
end

_product(::Type{ModularDouble{I}}) where {I} = ProductSector{Tuple{I,TimeReversed{I}}}

one_modular_double(::Type{I}) where {I<:Sector} = TensorKit.unit(ModularDouble{I})
doubleIsing(a::Symbol, b::Symbol)::ModularDouble{IsingAnyon} = ModularDouble(IsingAnyon(a), IsingAnyon(b))
function diagvalues(::Type{ModularDouble{I}})::Base.Generator where {I<:Sector}
    return (ModularDouble(a, a) for a in values(I))
end

function take_center(::Type{I})::Type{ModularDouble{I}} where {I<:Sector}
    TKS.ismodular(I) || throw(ArgumentError("no default center registered for $I"))
    return ModularDouble{I}
end

Base.isless(a::C,b::C) where {C<:ModularDouble} = isless(a.label,b.label)
Base.length(::TKS.SectorValues{C}) where {C<:ModularDouble} = length(values(_product(C)))
function Base.getindex(::TKS.SectorValues{C}, i::Int)::C where {C<:ModularDouble}
    return _from_product(values(_product(C))[i])
end
Base.iterate(v::TKS.SectorValues{C},i::Int=1) where {C<:ModularDouble} = i <= length(v) ? (v[i],i+1) : nothing
TKS.findindex(::TKS.SectorValues{C},a::C) where {C<:ModularDouble} = TKS.findindex(values(_product(C)),a.label)
function TensorKit.unit(::Type{ModularDouble{I}})::ModularDouble{I} where {I<:Sector}
    return ModularDouble(unit(I), unit(I))
end
function TensorKit.dual(a::C)::C where {C<:ModularDouble}
    return _from_product(dual(a.label))
end
function TensorKit.:⊗(a::C, b::C)::Base.Generator where {C<:ModularDouble}
    return (_from_product(c) for c in a.label ⊗ b.label)
end
for trait in (:FusionStyle,:BraidingStyle,:sectorscalartype,:fusionscalartype,:braidingscalartype)
    @eval TensorKit.$trait(::Type{C}) where {C<:ModularDouble} = TensorKit.$trait(_product(C))
end
for symbol in (:Nsymbol,:Rsymbol)
    @eval TensorKit.$symbol(a::C,b::C,c::C) where {C<:ModularDouble} = TensorKit.$symbol(a.label,b.label,c.label)
end
TensorKit.Fsymbol(a::C,b::C,c::C,d::C,e::C,f::C) where {C<:ModularDouble} =
    TensorKit.Fsymbol(a.label,b.label,c.label,d.label,e.label,f.label)
TensorKit.dim(a::ModularDouble) = dim(a.label)
TensorKit.twist(a::ModularDouble) = TensorKit.twist(a.label)

function _factors(a::ModularDouble{I}) where {I}
    l, r = a.label.sectors
    return Vect[I](l=>1) ⊗ Vect[I](r.a=>1)
end
function _boundary_labels(x::ModularDouble{I})::TKS.SectorSet{I} where {I<:Sector}
    left, right = x.label.sectors
    return TKS.SectorSet{I}(left ⊗ right.a)
end
function forget_Nsymbol(x::ModularDouble{I}, a::I)::Int where {I<:Sector}
    left, right = x.label.sectors
    return Int(Nsymbol(left, right.a, a))
end

# Slide the left factor over V and the time-reversed factor under V.
function HalfBraiding(a::ModularDouble{I}, V::ElementarySpace; T::Type{<:Number}=ComplexF64) where {I}
    sectortype(V) === I || throw(SectorMismatch("expected a boundary leg graded by $I"))
    factors = _factors(a)
    l, r = factors
    F = unitary(T, factors ← fuse(factors))
    @planar h[v x; y w] := F'[x; p q] * BraidingTensor(l,V)[v p; i k] *
        BraidingTensor(V,r)'[k q; j w] * F[i j; y]
    return h
end
function halfbraiding_symbol(a::ModularDouble{I}, i::I,r::I,s::I,j::I,μ::Int,ν::Int)::AbstractMatrix{<:Number} where {I}
    h = HalfBraiding(a,Vect[I](i=>1))
    f = TensorKit.FusionTree((i,r),j,(false,false),(),(μ,))
    g = TensorKit.FusionTree((s,i),j,(false,false),(),(ν,))
    return h[f,g][1,:,:,1]
end

# A coefficient-one splitting vertex in TensorKit's fixed fusion-tree basis.
function _split_vertex(a::I, b::I, c::I, channel::Int) where {I<:Sector}
    V = Vect[I]
    t = zeros(TensorKit.sectorscalartype(I), V(a=>1) ⊗ V(b=>1) ← V(c=>1))
    f = TensorKit.FusionTree((a,b),c,(false,false),(),(channel,))
    g = only(fusiontrees(domain(t),c))
    t[f,g] .= 1
    return t
end

function forgetvertexscalartype(::Type{ModularDouble{I}})::Type{<:Number} where {I<:Sector}
    return complex(float(TensorKit.sectorscalartype(I)))
end

function forget_vertex(a::C,b::C,c::C,r::I,s::I,t::I,λ::Int,μ::Int)::AbstractArray{<:Number,3} where {I,C<:ModularDouble{I}}
    al, ar = a.label.sectors
    bl, br = b.label.sectors
    cl, cr = c.label.sectors
    nright = TensorKit.Nsymbol(ar,br,cr)
    1 <= λ <= TensorKit.Nsymbol(a,b,c) || throw(ArgumentError("invalid bulk fusion channel $λ"))
    left, right = divrem(λ-1,nright)
    L = _split_vertex(al,bl,cl,left+1)
    R = _split_vertex(ar.a,br.a,cr.a,right+1)
    E = forgetvertexscalartype(C)
    A = unitary(E, _factors(a) ← forget(a))
    B = unitary(E, _factors(b) ← forget(b))
    Cmap = unitary(E, _factors(c) ← forget(c))
    cross = BraidingTensor(codomain(L)[2],codomain(R)[1])
    @planar image[x y; z] := A'[x; p q] * B'[y; u v] * cross[q u; k h] *
        L[p k; l] * R[h v; d] * Cmap[l d; z]
    f = TensorKit.FusionTree((r,s),t,(false,false),(),(μ,))
    g = only(fusiontrees(domain(image),t))
    return image[f,g]
end
