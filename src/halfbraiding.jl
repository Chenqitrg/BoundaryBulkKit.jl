function HalfBraiding_charge(a::QDℤ{N}, V::GradedSpace{ZNIrrep{N}, NTuple{N, Int64}}) where {N}
    fgt_a = forget_flux(a)
    W = ZNSpace{N}(fgt_a=>1)
    Ω = zeros(ComplexF64, V⊗W←W⊗V)
    for tree in fusiontrees(Ω)
        charge = tree[1].uncoupled[1]
        d = dim(V, charge)
        Ω[tree...] .= cis(2 * pi/N * a.flux * charge.n) * reshape(Matrix(I, d, d), (d, 1, 1, d))
    end
    return Ω
end

# Example:

# for a in SectorValues{QDℤ{2}}()
#     @show HalfBraiding_charge(a, Z2Space(0=>2,1=>2))
# end

# VecGIrr{ℤ₂}(ℤ₂(1))
# ℤ₂×ℤ₂
# QDAb{ℤ₂×ℤ₂}