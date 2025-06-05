

# function forget(a::ProductSector{T}) where {T<:Tuple{Vararg{QDℤ}}}
#     forgets = ()
#     for layer in a
#         forgets = (forgets..., forget(layer))
#     end
#     return fuse(⊠(forgets...))
# end

# function forget(a::ProductSector{Tuple{𝒞, TimeReversed{𝒞}}}) where {𝒞<:ModularSector}
#     aup, adown = a.sectors
#     return Vect[𝒞](aup=>1)⊗Vect[𝒞](adown.a=>1)
# end

# function _half_braiding_phase(bulk::ProductSector{T}, boundary::ProductSector{F}) where {T<:Tuple{Vararg{<:QDℤ}}, F<:Tuple{Vararg{<:ZNIrrep}}}
#     phase = 1.0 + 0 * im
#     for (anyon, charge) in zip(bulk, boundary)
#         N = typeof(charge).parameters[1]
#         phase *= cis(2 * pi/N * anyon.flux * charge.n)
#     end
#     return phase
# end

# function HalfBraiding(a::ProductSector{T}, V::GradedSpace{ProductSector{F}, G}) where {T<:Tuple{Vararg{<:QDℤ}}, F<:Tuple{Vararg{<:ZNIrrep}}, G<:Tuple{Vararg{Int}}}
#     if length(a.sectors) != length(F.parameters)
#         throw(ArgumentError("The layer of object $a does not match that of the graded space $V"))
#     end

#     fgt_a = forget(a)

#     if sectortype(fgt_a) != sectortype(V)
#         throw(ArgumentError("The sector of forget($a) does not match that of the graded space $V"))
#     end

#     Ω = zeros(ComplexF64, V⊗fgt_a←fgt_a⊗V)

#     for tree in fusiontrees(Ω)
#         charge = tree[1].uncoupled[1]
#         d = dim(V, charge)
#         Ω[tree...] .= _half_braiding_phase(a, charge) * reshape(Matrix(I, d, d), (d, 1, 1, d))
#     end
#     return Ω
# end

# function HalfBraiding(a::ProductSector{Tuple{𝒞, TimeReversed{𝒞}}}, V::GradedSpace{𝒞, T}) where {𝒞<:ModularSector, T<:Tuple{Vararg{Int}}}
#     WW = forget(a)
#     Wup, Wdown = WW.spaces
#     W = fuse(WW)
#     WWTW = unitary(WW ← W)
#     @planar Ω[vu wu; wd vd] := WWTW'[wu; newup newdown] * BraidingTensor(Wup, V)[vu newup; up vmd] * BraidingTensor(V, Wdown)'[vmd newdown; down vd] * WWTW[up down; wd]
#     return Ω
# end


# Example:

# for a in SectorValues{QDℤ{2}}()
#     @show HalfBraiding_charge(a, Z2Space(0=>2,1=>2))
# end

# VecGIrr{ℤ₂}(ℤ₂(1))
# ℤ₂×ℤ₂
# QDAb{ℤ₂×ℤ₂}
