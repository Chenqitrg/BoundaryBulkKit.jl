include(joinpath(@__DIR__, "support", "setup.jl"))

let
    # In Rep(A₄), the triplet appears twice in triplet ⊗ triplet.
    r = A4Irrep(3)
    a = ElectricA4(r)
    V, W = Vect[A4Irrep](r => 1), Vect[ElectricA4](a => 1)
    U = zeros(Int, CentralTensorMapSpace(V ⊗ V, W ⊗ W))
    block(U, a) .= [11 12; 21 22]

    for (fout, v, fin) in fusiontrees(U)
        v.bulk == a || continue
        row, col = only(fout.vertices), only(fin.vertices)
        value = only(U[fout, v, fin])
        @show row col value
        @assert value == 10 * row + col
    end
end
