# Only the forgetful decomposition is needed to test junctions with multiplicity.
struct JunctionFixture <: DrinfeldCenter{Z2Irrep} end
BoundaryBulkKit.forget(::JunctionFixture) = Vect[Z2Irrep](Z2Irrep(0) => 1, Z2Irrep(1) => 2)

@testset "Bulk–boundary junction basis" begin
    a = JunctionFixture()
    even, odd = Z2Irrep(0), Z2Irrep(1)
    v₁, v₂ = junctions(a, odd)
    @test (v₂.bulk, v₂.boundary, v₂.copyindex) == (a, odd, 2)
    @test Set(junctions(a)) == Set((BoundaryBulkJunction(a, even, 1), v₁, v₂))
    @test isempty(junctions(Z2Center(true, false), even))
    @test_throws BoundsError BoundaryBulkJunction(a, odd, 0)
    @test_throws BoundsError BoundaryBulkJunction(a, odd, 3)
    @test_throws BoundsError BoundaryBulkJunction(Z2Center(true, false), even, 1)
    @test_throws MethodError BoundaryBulkJunction(a, Z3Irrep(1), 1)

    P = TensorMap(v₂)
    @test space(P) == (Vect[Z2Irrep](odd => 1) ← forget(a))
    @test block(P, odd) == [0 1]
    @test sum(v -> TensorMap(v)' * TensorMap(v), junctions(a)) ≈ id(forget(a))
end
