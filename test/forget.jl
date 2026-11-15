@testset "Forgetful functor" begin
    C = Z3QuantumDouble
    x, y = C(1, 0), C(1, 1)
    V = Vect[C](x => 2, y => 1)
    @test forget_Nsymbol(x, Z3Irrep(1)) === 1
    @test forget_Nsymbol(x, Z3Irrep(0)) === 0
    A = zeros(ComplexF64, V ← V)
    block(A, x) .= [1 2im; 3 4]
    block(A, y) .= 5
    @test block(@inferred(forget(A)), Z3Irrep(1)) == [1 2im 0; 3 4 0; 0 0 5]
    @test forget(A ⊗ id(V')) ≈ forget(A) ⊗ id(forget(V'))
    @test @inferred(forget(id(Vect[C]()))) == id(forget(Vect[C]()))

    J = Vect[JunctionFixture](JunctionFixture() => 2)
    D = zeros(ComplexF64, J ← J)
    block(D, JunctionFixture()) .= [1 2; 3 4]
    @test block(forget(D), Z2Irrep(1)) == kron([1 2; 3 4], [1 0; 0 1])

    for C in (take_center(Trivial), ModularDouble{FibonacciAnyon})
        a = last(collect(values(C)))
        X = Vect[C](a => 1)
        P = X ⊗ X ⊗ X
        W = forget(X) ⊗ forget(X) ⊗ forget(X)
        A = randn(ComplexF64, P ← X ⊗ X)
        B = randn(ComplexF64, X ⊗ X ← X)
        @test forget(id(P)) ≈ id(W)
        @test @inferred(forget(A * B)) ≈ forget(A) * forget(B)
        @test @inferred(forget(A')) ≈ forget(A)'
        @test forget(B ⊗ id(X)) ≈ forget(B) ⊗ id(forget(X))
        state = randn(ComplexF64, P ← one(P))
        @test forget(state') * forget(state) ≈ forget(state' * state)
    end
end

@testset "Forgetful scalar type and inference" begin
    for (C, expected) in ((Z3QuantumDouble, Float32), (Z2Center, Float32),
                          (take_center(Trivial), Float32), (ModularDouble{FibonacciAnyon}, ComplexF64))
        V = Vect[C](unit(C) => 1)
        for P in (V ⊗ V, V ⊗ Vect[C]())
            image = @inferred forget(zeros(Float32, P ← V))
            @test scalartype(image) === expected
        end
        scalar = zeros(Float32, one(ProductSpace(V)) ← one(ProductSpace(V)))
        @test scalartype(@inferred forget(scalar)) === Float32
    end
end
