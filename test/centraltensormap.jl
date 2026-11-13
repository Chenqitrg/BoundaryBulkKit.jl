@testset "CentralTensorMap sectors and shared storage" begin
    e, m = Z2Center(true, false), Z2Center(false, true)
    vacuum, odd = Z2Irrep(0), Z2Irrep(1)
    V = Vect[Z2Irrep](vacuum => 1, odd => 2)
    W = Vect[Z2Center](e => 2, m => 1)
    S = CentralTensorMapSpace(V ⊗ V, W ⊗ W)
    @test S isa CentralTensorMapSpace{Z2Irrep,Z2Center}
    @test_throws TypeError CentralTensorMapSpace{Z3Irrep,Z2Center}
    # The unit and em channels each have five or four boundary rows.
    @test dim(S) == 5 * (2 * 2 + 1) + 4 * (2 + 2)
    data = collect(1:dim(S))
    U = CentralTensorMap(data, S)
    @test U.data === data
    @test space(U) === S
    @test (codomain(U), domain(U)) == (V ⊗ V, W ⊗ W)
    @test (numout(U), numin(U), scalartype(U)) == (2, 2, Int)
    @test size(block(U, unit(Z2Center))) == (5, 5)
    @test size(block(U, only(e ⊗ m))) == (4, 4)
    @test size(block(U, m)) == (5, 0)
    @test Set(first.(collect(blocks(U)))) == Set(blocksectors(U))
    @test_throws DimensionMismatch CentralTensorMap(Int[], S)
    @test_throws TensorKit.SectorMismatch CentralTensorMapSpace(W, V)
    @test_throws TensorKit.SectorMismatch CentralTensorMapSpace(Vect[Z3Irrep](Z3Irrep(0) => 1), W)

    @test length(subblocks(U)) == length(fusiontrees(U)) == 8
    seen = Int[]
    for (f, v, g) in fusiontrees(U)
        A = U[f, v, g]
        @test size(A) == (TensorKit.dims(codomain(U), f.uncoupled)..., TensorKit.dims(domain(U), g.uncoupled)...)
        append!(seen, vec(A))
    end
    @test sort(seen) == data

    a = unit(Z2Center)
    f = only(f for f in fusiontrees(codomain(U), vacuum) if f.uncoupled == (odd, odd))
    g = only(g for g in fusiontrees(domain(U), a) if g.uncoupled == (e, e))
    v = only(junctions(a))
    A = U[f, v, g]
    @test size(A) == (2, 2, 2, 2)
    # The (m,m) tree occupies column 1. The four (e,e) columns follow it.
    @test A[2, 1, 2, 1] == block(U, a)[3, 3]
    A[2, 1, 2, 1] = -8
    @test data[13] == block(U, a)[3, 3] == -8
    block(U, a)[3, 3] = -12
    @test A[2, 1, 2, 1] == -12
    @test_throws TensorKit.SectorMismatch U[f, only(junctions(m)), g]
    @test_throws TensorKit.SectorMismatch block(U, odd)
end

@testset "CentralTensorMap one bulk leg and empty spaces" begin
    a = Z2Center(true, false)
    r = Z2Irrep(1)
    V = Vect[Z2Irrep](r => 2)
    W = Vect[Z2Center](a => 3)
    U = CentralTensorMap{ComplexF64}(undef, V, W)
    block(U, a) .= reshape(1:6, 2, 3)
    @test U[only(fusiontrees(U))...] == reshape(1:6, 2, 3)
    empty = CentralTensorMap{Float64}(undef, Vect[Z2Irrep](), W)
    @test isempty(empty.data) && isempty(blocksectors(empty)) && isempty(fusiontrees(empty))
    @test size(block(empty, a)) == (0, 3)
end

# Supply the single-leg sector interface for the multiplicity fixture.
TensorKit.unit(::Type{JunctionFixture}) = JunctionFixture()
TensorKit.dual(a::JunctionFixture) = a
TensorKit.FusionStyle(::Type{JunctionFixture}) = UniqueFusion()
Base.length(::typeof(values(JunctionFixture))) = 1
Base.getindex(::typeof(values(JunctionFixture)), i::Int) = (JunctionFixture(),)[i]
Base.iterate(::typeof(values(JunctionFixture)), i::Int=1) = iterate((JunctionFixture(),), i)

@testset "Junction copies occupy independent rows in TensorKit state order" begin
    a = JunctionFixture()
    even, odd = Z2Irrep(0), Z2Irrep(1)
    V = Vect[Z2Irrep](even => 1, odd => 2)
    W = Vect[JunctionFixture](a => 2)
    ψ = zeros(Float64, V ⊗ V ← forget(a))
    ψ.data .= eachindex(ψ.data)
    U = CentralTensorMap{Float64}(undef, codomain(ψ), W)
    B = block(U, a)
    @test size(B) == (length(ψ.data), 2)
    B[:, 1] .= ψ.data
    B[:, 2] .= -ψ.data
    for (f, v, g) in fusiontrees(U)
        _, h = only((p, q) for (p, q) in fusiontrees(ψ) if p == f)
        expected = view(ψ[f, h], :, :, v.copyindex)
        A = U[f, v, g]
        @test A[:, :, 1] == expected
        @test A[:, :, 2] == -expected
    end
end

@testset "View access through TensorKit entry points" begin
    a = Z2Center(true, false)
    V = Vect[Z2Irrep](Z2Irrep(1) => 2)
    W = Vect[Z2Center](a => 2)
    U = CentralTensorMap(zeros(4), V', W')
    k = only(fusiontrees(U))
    U[k...] = [1 2; 3 4]
    @test eltype(U) === scalartype(typeof(U)) === Float64
    @test TensorKit.storagetype(U) === Vector{Float64}
    @test k[1].isdual == k[3].isdual == (true,)
    @test blocks(U)[a] == [1 2; 3 4]
    subblocks(U)[k][1, 2] = 9
    @test block(U, a)[1, 2] == 9
    @test only(subblocks(U)).first == k

    S = CentralTensorMapSpace(one(V), one(W))
    scalar = zeros(Int, S)
    k₀ = only(fusiontrees(scalar))
    scalar[k₀...][] = 7
    @test scalar.data == [7]
    @test size(scalar[k₀...]) == ()
end
