@testset "Group and modular center adapters" begin
    C = take_center(Z3Irrep)
    @test C <: QuantumDouble{ℤ{3}} <: DrinfeldCenter
    groupof(::QuantumDouble{G}) where {G} = G
    @test groupof(C(1, 0)) === ℤ{3}
    @test C <: DrinfeldCenter{Z3Irrep}
    @test C === 𝒟Z3 === 𝒟ZN{3}
    @test !(C <: QuantumDouble{ℤ{2}})
    @test base_sector(C) === Z3Irrep
    @test all(i -> TensorKit.findindex(values(C),values(C)[i]) == i, 1:9)
    @test only(halfbraiding_symbol(C(1,1),Z3Irrep(1),Z3Irrep(1),Z3Irrep(1),Z3Irrep(2),1,1)) ≈ cispi(2/3)
    F = take_center(FibonacciAnyon)
    @test F === ModularDouble{FibonacciAnyon}
    @test !(F <: QuantumDouble)
    @test first(diagvalues(F)) == one_modular_double(FibonacciAnyon)
    @test doubleIsing(:σ, :σ) in diagvalues(take_center(IsingAnyon))
    a = last(collect(values(F)))
    V = Vect[FibonacciAnyon](:τ=>1)
    @test HalfBraiding(a,V) ≈ HalfBraiding(Vect[F](a=>1),V)
    for c in a ⊗ a
        p = zeros(ComplexF64,forget(a) ⊗ forget(a) ← forget(c))
        for (f,g) in fusiontrees(p)
            p[f,g] .= forget_vertex(a,a,c,f.uncoupled...,g.coupled,1,only(f.vertices))
        end
        @test p' * p ≈ id(forget(c))
    end
end
