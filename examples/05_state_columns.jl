include(joinpath(@__DIR__, "support", "setup.jl"))
using LinearAlgebra: Symmetric, eigen

let
    a = Z2QuantumDouble(1, 0)
    r = Z2Irrep(1)
    V = Vect[Z2Irrep](r => 2)
    T = zeros(Float64, V ← V)
    block(T, r) .= [2 1; 1 2]

    # Diagonalize the small boundary block, then store the two right states.
    D, vectors = eigen(Symmetric(block(T, r)))
    U = CentralTensorMap{Float64}(undef, V, Vect[Z2QuantumDouble](a => 2))
    B = block(U, a)
    for j in eachindex(D)
        ψ = zeros(Float64, V ← forget(a))
        block(ψ, r)[:, 1] .= vectors[:, j]
        B[:, j] .= ψ.data
    end

    # Read the second state back as a TensorKit map and verify its eigenvalue.
    ψ₂ = TensorMap{Float64}(B[:, 2], V ← forget(a))
    @show D block(ψ₂, r)
    @assert T * ψ₂ ≈ D[2] * ψ₂
end
