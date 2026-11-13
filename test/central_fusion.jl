@testset "Independent boundary and bulk fusion paths" begin
    I = FibonacciAnyon
    C = take_center(I)
    τ = I(:τ)
    a = C(τ, TimeReversed{I}(unit(I)))
    V, W = Vect[I](τ => 1), Vect[C](a => 1)
    U = zeros(Int, CentralTensorMapSpace(V ⊗ V ⊗ V, W ⊗ W ⊗ W))
    B = block(U, a)
    @test size(B) == (2, 2)
    B .= [11 12; 21 22]
    rows = collect(fusiontrees(codomain(U), τ))
    cols = collect(fusiontrees(domain(U), a))
    for (row, f) in enumerate(rows), (col, g) in enumerate(cols)
        @test only(U[f, only(junctions(a)), g]) == B[row, col]
    end
end
