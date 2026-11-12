include(joinpath(@__DIR__, "support", "setup.jl"))

let
    a = Z2Center(true, false)  # Electric charge.
    r = Z2Irrep(1)
    V = Vect[Z2Irrep](r => 2)
    W = Vect[Z2Center](a => 3)

    data = collect(1:6)
    U = CentralTensorMap(data, V, W)
    B = block(U, a)
    @show data B
    @assert B == [1 3 5; 2 4 6]

    B[2, 3] = 99
    @show data
    @assert data[6] == 99
end
