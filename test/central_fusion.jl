@testset "Independent boundary and bulk fusion channels" begin
    r = A4Irrep(3)
    a = ElectricA4(r)
    V, W = Vect[A4Irrep](r => 1), Vect[ElectricA4](a => 1)
    U = CentralTensorMap{Int}(undef, V ⊗ V, W ⊗ W)
    B = block(U, a)
    @test size(B) == (2, 2)
    B .= [1 3; 2 4]
    for (f, v, g) in fusiontrees(U)
        v.bulk == a || continue
        @test only(U[f, v, g]) == B[only(f.vertices), only(g.vertices)]
    end
end
