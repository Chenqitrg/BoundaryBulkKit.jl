include(joinpath(@__DIR__, "support", "setup.jl"))
if !isdefined(@__MODULE__, :ExampleCenters)
    include(joinpath(@__DIR__, "support", "known_centers.jl"))
end
using .ExampleCenters: Z2Center

let
    # Compare the independently implemented Z₂ center with the built-in model.
    @assert base_sector(Z2Center) === Z2Irrep
    @assert TensorKit.ismodular(Z2Center)
    V = Vect[Z2Irrep](0 => 1, 1 => 2)
    for a in values(Z2Center)
        b = Z2QuantumDouble(a.charge, a.flux)
        @assert forget(a) == forget(b)
        @assert HalfBraiding(a, V) ≈ HalfBraiding(b, V)
    end
    e, m = Z2Center(true, false), Z2Center(false, true)
    c = only(e ⊗ m)
    r, s, t = (only(sectors(forget(x))) for x in (e, m, c))
    @assert forget_vertex(e, m, c, r, s, t, 1, 1) == ones(Int, 1, 1, 1)
end
