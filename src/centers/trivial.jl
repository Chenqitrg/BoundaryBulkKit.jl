"""The unique simple object of the center of plain vector spaces."""
struct TrivialCenter <: DrinfeldCenter{Trivial} end
take_center(::Type{Trivial}) = TrivialCenter
TensorKit.unit(::Type{TrivialCenter}) = TrivialCenter()
TensorKit.dual(a::TrivialCenter) = a
TensorKit.FusionStyle(::Type{TrivialCenter}) = TensorKit.UniqueFusion()
TensorKit.BraidingStyle(::Type{TrivialCenter}) = TensorKit.Bosonic()
TensorKit.:⊗(a::TrivialCenter, ::TrivialCenter) = (a,)
TensorKit.Nsymbol(::TrivialCenter,::TrivialCenter,::TrivialCenter) = true
TensorKit.Fsymbol(::Vararg{TrivialCenter,6}) = 1
TensorKit.Rsymbol(::Vararg{TrivialCenter,3}) = 1
Base.isless(::TrivialCenter,::TrivialCenter) = false
Base.length(::TKS.SectorValues{TrivialCenter}) = 1
Base.getindex(::TKS.SectorValues{TrivialCenter},i::Int) = (TrivialCenter(),)[i]
Base.iterate(::TKS.SectorValues{TrivialCenter}, i::Int=1) = iterate((TrivialCenter(),), i)
TKS.findindex(::TKS.SectorValues{TrivialCenter},::TrivialCenter) = 1
forget(::TrivialCenter) = TensorKit.ComplexSpace(1)
halfbraiding_symbol(::TrivialCenter,::Trivial,::Trivial,::Trivial,::Trivial,::Int,::Int) = ones(Int,1,1)
forget_vertex(::TrivialCenter,::TrivialCenter,::TrivialCenter,::Trivial,::Trivial,::Trivial,::Int,::Int) = ones(Int,1,1,1)
