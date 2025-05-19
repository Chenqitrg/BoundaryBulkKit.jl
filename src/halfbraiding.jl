forget_flux(a::QDℤ{N}) where {N} = ZNSpace{N}(a.charge=>1)
foget_charge(a::QDℤ{N}) where {N} = ZNSpace{N}(a.flux=>1)

function unfolding(a::ProductSector{Tuple{𝒞, TimeReversed{𝒞}}}) where {𝒞<:ModularSector}
    aup, adown = a.sectors
    return Vect[𝒞](aup=>1), Vect[𝒞](adown.a=>1)
end

function HalfBraiding(a::QDℤ{N}, V::GradedSpace{ZNIrrep{N}, NTuple{N, Int64}}) where {N}
    fgt_a = forget_flux(a)
    Ω = zeros(ComplexF64, V⊗fgt_a←fgt_a⊗V)
    for tree in fusiontrees(Ω)
        charge = tree[1].uncoupled[1]
        d = dim(V, charge)
        Ω[tree...] .= cis(2 * pi/N * a.flux * charge.n) * reshape(Matrix(I, d, d), (d, 1, 1, d))
    end
    return Ω
end

function HalfBraiding(a::ProductSector{Tuple{𝒞, TimeReversed{𝒞}}}, V::GradedSpace{𝒞, T}) where {𝒞<:ModularSector, T<:Tuple{Vararg{Int}}}
    Wup, Wdown = unfolding(a)
    W = fuse(Wup ⊗ Wdown)
    WWTW = unitary(Wup ⊗ Wdown ← W)
    @planar Ω[vu wu; wd vd] := WWTW'[wu; newup newdown] * BraidingTensor(Wup, V)[vu newup; up vmd] * BraidingTensor(V, Wdown)'[vmd newdown; down vd] * WWTW[up down; wd]
    return Ω
end
# Example:

# for a in SectorValues{QDℤ{2}}()
#     @show HalfBraiding_charge(a, Z2Space(0=>2,1=>2))
# end

# VecGIrr{ℤ₂}(ℤ₂(1))
# ℤ₂×ℤ₂
# QDAb{ℤ₂×ℤ₂}