include(joinpath(@__DIR__, "support", "setup.jl"))

let
    C = take_center(Trivial)
    a = unit(C)
    X = forget(a)
    @assert X == ℂ^1
    V = ℂ^2
    @assert HalfBraiding(a, V) ≈ TensorMap(BraidingTensor(X, V))
    U = CentralTensorMap(collect(1:4), V, Vect[C](a => 2))
    @assert block(U, a) == [1 3; 2 4]
end
