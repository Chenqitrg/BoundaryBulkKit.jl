using TensorKit
using TensorKitSectors
using LinearAlgebra

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

HalfBraiding_charge(QDℤ{2}(1,1), Z2Space(0=>2,1=>2))
# V = Z2Space(0=>2,1=>2)
# dim(V, Z2Irrep(0))