@testset "Half-braiding a graded center leg" begin
    even, odd = Z2Irrep(0), Z2Irrep(1)
    e, m, em = Z2Center(true, false), Z2Center(false, true), Z2Center(true, true)
    A = Vect[Z2Center](m => 1, e => 2, em => 1)
    V = Vect[Z2Irrep](even => 1, odd => 2)
    X = Vect[Z2Irrep](even => 1, odd => 3)
    @test forget(A) == X
    for center in (A, A'), boundary in (V, V')
        Γ = HalfBraiding(center, boundary)
        W = forget(center)
        @test W == (isdual(center) ? X' : X)
        @test space(Γ) == (boundary ⊗ W ← W ⊗ boundary)
        @test Γ' * Γ ≈ id(domain(Γ))
        for (f, g) in fusiontrees(Γ)
            i, r = f.uncoupled
            s, k = g.uncoupled
            B = Γ[f, g]
            for p in CartesianIndices(B)
                b, α, β, d = Tuple(p)
                flux = r == even || α == 3
                expected = (i == k && r == s && b == d && α == β) * (flux && i == odd ? -1 : 1)
                @test B[p] == expected
            end
        end
    end
    @test HalfBraiding(e, V) == HalfBraiding(Vect[Z2Center](e => 1), V)
    @test isempty(HalfBraiding(Vect[Z2Center](), V).data)
    @test isempty(HalfBraiding(A, Vect[Z2Irrep]()).data)
    @test_throws SectorMismatch HalfBraiding(A, Vect[Z3Irrep](Z3Irrep(1) => 1))
end

# The electric embedding uses the ordinary Rep(A₄) crossing.
function BoundaryBulkKit.halfbraiding_symbol(a::ElectricA4, i::A4Irrep, r::A4Irrep,
        s::A4Irrep, j::A4Irrep, μ::Int, ν::Int)
    @assert r == s == a.rep
    return fill(Rsymbol(s, i, j)[ν, μ], 1, 1)
end

@testset "Dual non-self-dual charges and fusion multiplicity" begin
    A = Vect[ElectricA4](ElectricA4(A4Irrep(1)) => 2, ElectricA4(A4Irrep(3)) => 1)
    V = Vect[A4Irrep](A4Irrep(2) => 1, A4Irrep(3) => 2)
    @test forget(A') == forget(A)'
    for center in (A, A'), boundary in (V, V')
        Γ = HalfBraiding(center, boundary)
        @test Γ ≈ TensorMap(BraidingTensor(forget(center), boundary))
    end
end

# Use an involution on a pair of junction copies to test complex matrix
# elements independently of the bulk-copy and boundary-copy identities.
function BoundaryBulkKit.halfbraiding_symbol(::JunctionFixture, i::Z2Irrep, r::Z2Irrep,
        s::Z2Irrep, j::Z2Irrep, μ::Int, ν::Int)
    @assert r == s && j == only(i ⊗ r) && μ == ν == 1
    r == Z2Irrep(0) && return ones(1, 1)
    return i == Z2Irrep(1) ? [0 im; -im 0] : [1 0; 0 1]
end

@testset "Complex junction matrix inside each copy" begin
    a, odd = JunctionFixture(), Z2Irrep(1)
    A = Vect[JunctionFixture](a => 2)
    V = Vect[Z2Irrep](odd => 3)
    Γ = HalfBraiding(A, V; T=ComplexF32)
    @test scalartype(Γ) === ComplexF32
    f, g = only((f, g) for (f, g) in fusiontrees(Γ) if f.uncoupled == (odd, odd))
    B = Γ[f, g]
    @test size(B) == (3, 4, 4, 3)
    @test B[2, :, :, 2] == [0 im 0 0; -im 0 0 0; 0 0 0 im; 0 0 -im 0]
    @test iszero(B[1, :, :, 2])
end
