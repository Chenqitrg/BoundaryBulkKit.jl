"""
    struct QDℤ{N} <: Sector
        charge::Int
        flux::Int
    end
    
The quantum double of Rep(ℤₙ) or Vec(ℤₙ).
"""
struct QDℤ{N} <: Sector
    charge::Int
    flux::Int
    function QDℤ{N}(e::Int, m::Int) where {N}
        new{N}(mod(e, N), mod(m, N))
    end
end

ℨ(::Type{ZNIrrep{N}}) where {N} = QDℤ{N}

FusionStyle(::Type{QDℤ{N}}) where {N} = SimpleFusion()
BraidingStyle(::Type{QDℤ{N}}) where {N} = Anyonic()

Nsymbol(a::QDℤ{N}, b::QDℤ{N}, c::QDℤ{N}) where {N} =
    (mod(a.charge+b.charge, N)==c.charge) && (mod(a.flux+b.flux, N)==c.flux)
Fsymbol(a::QDℤ{N}, b::QDℤ{N}, c::QDℤ{N}, d::QDℤ{N}, e::QDℤ{N}, f::QDℤ{N}) where {N} =
    Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d)
Base.one(::Type{QDℤ{N}}) where {N} = QDℤ{N}(0, 0)
Base.conj(c::QDℤ{N}) where {N} = QDℤ{N}(-c.charge, -c.flux)
⊗(c1::QDℤ{N}, c2::QDℤ{N}) where {N} = (QDℤ{N}(c1.charge+c2.charge, c1.flux+c2.flux),)
Base.IteratorSize(::Type{SectorValues{QDℤ{N}}}) where {N} = HasLength()
Base.length(::SectorValues{QDℤ{N}}) where {N} = N^2
function Base.getindex(::SectorValues{QDℤ{N}}, i::Int) where {N}
    e = (i-1)÷N
    m = (i-1)%N
    return QDℤ{N}(e, m)
end
Base.iterate(::SectorValues{QDℤ{N}}, i::Int = 1) where {N} =
    i == N^2+1 ? nothing : (QDℤ{N}((i-1)÷N, (i-1)%N), i + 1)
findindex(::SectorValues{QDℤ{N}}, a::QDℤ{N}) where {N} = a.charge*N + a.flux + 1
function Base.isless(c1::QDℤ{N}, c2::QDℤ{N}) where {N}
    if c1.charge < c2.charge
        return true
    elseif (c1.charge == c2.charge) && c1.flux < c2.flux
        return true
    else
        return false
    end
end
function Rsymbol(a::QDℤ{N}, b::QDℤ{N}, c::QDℤ{N}) where {N}
    R = Nsymbol(a, b, c) * exp(2 * pi * im/N * a.flux * b.charge)
    return R
end

forget(a::QDℤ{N}) where {N} = ZNSpace{N}(a.charge=>1)
function HalfBraiding(a::QDℤ{N}, V::GradedSpace{ZNIrrep{N},NTuple{N,Int64}}) where {N}
    fgt_a = forget(a)
    Ω = zeros(ComplexF64, V⊗fgt_a←fgt_a⊗V)
    for tree in fusiontrees(Ω)
        charge = tree[1].uncoupled[1]
        d = dim(V, charge)
        Ω[tree...] .=
            cis(2 * pi/N * a.flux * charge.n) * reshape(Matrix(I, d, d), (d, 1, 1, d))
    end
    return Ω
end