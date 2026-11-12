include(joinpath(@__DIR__, "support", "setup.jl"))

let
    e, m = Z2Center(true, false), Z2Center(false, true)
    em = only(e ⊗ m)
    V = forget(em)
    W₁ = Vect[Z2Center](e => 2)
    W₂ = Vect[Z2Center](m => 3)
    U = CentralTensorMap(collect(1:6), V, W₁ ⊗ W₂)

    fout, v, fin = only(fusiontrees(U))
    A = U[fout, v, fin]
    @show fin.uncoupled fin.coupled block(U, em) size(A) A[1, :, :]
    @assert fin.uncoupled == (e, m) && fin.coupled == em
    @assert A[1, :, :] == [1 3 5; 2 4 6]
end
