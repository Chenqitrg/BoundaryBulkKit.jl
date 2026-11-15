include(joinpath(@__DIR__, "support", "setup.jl"))

let
    e  = Z2QuantumDouble(1, 0)
    em = Z2QuantumDouble(1, 1)
    @assert forget(e) == forget(em)  # Both have odd boundary charge.

    V = Vect[Z2Irrep](Z2Irrep(1) => 2)
    W = Vect[Z2QuantumDouble](e => 1, em => 2)
    U = zeros(Int, CentralTensorMapSpace(V, W))
    block(U, e)[:, 1] .= [1, 2]
    block(U, em) .= [3 5; 4 6]

    for (a, B) in blocks(U)
        @show a B
    end
    @show U.data
    @assert U.data == [1, 2, 3, 4, 5, 6]
end
