include(joinpath(@__DIR__, "support", "setup.jl"))

let
    I = FibonacciAnyon
    C = take_center(I)
    τ = I(:τ)
    a = C(τ, TimeReversed{I}(τ))
    X = forget(a)
    @assert X == Vect[I](unit(I) => 1, τ => 1)
    @assert length(collect(junctions(a))) == 2

    V = Vect[I](τ => 1)
    Γ = HalfBraiding(a, V)
    @assert space(Γ) == (V ⊗ X ← X ⊗ V)
    @assert Γ' * Γ ≈ id(domain(Γ))

    # Reconstruct the image of one bulk splitting vertex in boundary tree order.
    c = unit(C)
    p = zeros(ComplexF64, X ⊗ X ← forget(c))
    for (f, g) in fusiontrees(p)
        r, s = f.uncoupled
        p[f,g] .= forget_vertex(a, a, c, r, s, g.coupled, 1, only(f.vertices))
    end
    @assert p' * p ≈ id(forget(c))
end
