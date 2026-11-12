include(joinpath(@__DIR__, "support", "setup.jl"))

let
    e, em = Z2Center(true, false), Z2Center(true, true)
    even, odd = Z2Irrep(0), Z2Irrep(1)

    # Each object has one odd junction. Crossing an odd line detects the flux.
    for a in (e, em)
        v = only(junctions(a, odd))
        P = TensorMap(v)
        Ω = halfbraiding_symbol(a, odd, odd, odd, even, 1, 1)
        @show a v.copyindex block(P, odd) Ω
        @assert only(block(P, odd)) == 1
        @assert only(Ω) == (a == e ? 1 : -1)
    end

    # Cross a direct sum of both bulk charges with two copies of the odd leg.
    A = Vect[Z2Center](e => 2, em => 1)
    V = Vect[Z2Irrep](odd => 2)
    Γ = HalfBraiding(A, V)
    fout, fin = only(fusiontrees(Γ))
    @show space(Γ) Γ[fout, fin][1, :, :, 1]
    @assert Γ[fout, fin][1, :, :, 1] == [1 0 0; 0 1 0; 0 0 -1]
    @assert Γ' * Γ ≈ id(forget(A) ⊗ V)
end
