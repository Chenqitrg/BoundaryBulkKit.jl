include(joinpath(@__DIR__, "support", "setup.jl"))

let
    a = unit(Z2QuantumDouble)
    even, odd = Z2Irrep(0), Z2Irrep(1)
    V = Vect[Z2Irrep](even => 1, odd => 1)
    W = Vect[Z2QuantumDouble](a => 2)
    U = CentralTensorMap(collect(1:4), V ⊗ V, W)

    # The (even, even) and (odd, odd) output trees give one row each.
    for (fout, v, fin) in fusiontrees(U)
        A = U[fout, v, fin]
        @show fout.uncoupled v.boundary fin.uncoupled size(A) vec(A)
        if fout.uncoupled == (odd, odd)
            @assert vec(A) == [2, 4]
            A[1, 1, 2] = 40
        end
    end
    @show block(U, a)
    @assert block(U, a) == [1 3; 2 40]
end
