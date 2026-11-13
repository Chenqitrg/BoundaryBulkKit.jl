# BoundaryBulkKit.jl

BoundaryBulkKit uses TensorKit to represent simple objects of a Drinfeld center,
their forgetful images on the boundary, and tensor maps with both bulk and
boundary legs.

## Center families

| Boundary sector type `I` | `take_center(I)` | Worked example |
| --- | --- | --- |
| `Trivial` | The center with one simple object. Its forgetful image is `ℂ^1`. | [Trivial center and ordinary matrix](examples/09_trivial_center.jl) |
| `ZNIrrep{N}` | `ZNQuantumDouble{N}`, with labels `(charge, flux)`. Forgetting retains the charge. Crossing a boundary charge detects the flux. | [Junction projections and half-braiding](examples/06_junction_and_crossing.jl) |
| A modular sector such as `FibonacciAnyon` or `IsingAnyon` | `ModularDouble{I}`, with labels `(left, time-reversed right)`. A simple bulk object can decompose into several boundary sectors. | [Fibonacci decomposition, crossing and splitting vertex](examples/07_modular_double.jl) |

## Implementing another center

Choose a concrete `C <: DrinfeldCenter{I}` for the full center. Supply its
ordinary TensorKit sector operations (simple labels, unit, dual, fusion,
ordering, F and R symbols). BoundaryBulkKit supplies `base_sector(C) === I`
and the TensorKitSectors trait `ismodular(C) == true`. To make it the default
for a new boundary type, define `BoundaryBulkKit.take_center(::Type{I}) = C`.

The three additional methods use the same one-based junction basis. Lowercase
`a, b, c` below have type `C`; `i, j, r, s, t` have type `I`.
Write `m(a,r) = dim(forget(a),r)` for the number of copies.

| Method | Result | Meaning of one entry |
| --- | --- | --- |
| `forget(a)` | An `ElementarySpace` with sector type `I` | Its multiplicities specify the junctions from `a` to each boundary sector. |
| `halfbraiding_symbol(a,i,r,s,j,μ,ν)` | Numeric matrix of size `(m(a,r), m(a,s))` | Entry `[p,q]` maps copy `q` of `s` to copy `p` of `r` while crossing `i`. Channels `μ,ν` select the output and input trees of total charge `j`. |
| `forget_vertex(a,b,c,r,s,t,λ,μ)` | Numeric rank-three array of size `(m(a,r), m(b,s), m(c,t))` | Entry `[p,q,k]` maps copy `k` of `t` to copies `p,q` of `r,s`, for bulk channel `λ` and boundary channel `μ`. |

Read the [custom Z₂ implementation](examples/support/known_centers.jl) and run
[its comparison with the built-in quantum double](examples/08_custom_center.jl).
The two implementations can coexist; this example leaves the default for
`Z2Irrep` with the built-in center.

## Working with central tensor maps

Run any example from the repository root, for instance
`julia --project=. examples/01_one_block.jl`. Each script constructs its inputs
and checks its result.

| Task | Examples |
| --- | --- |
| Build and edit sector blocks | [One block and its storage](examples/01_one_block.jl) · [Two bulk sectors with the same boundary charge](examples/02_two_bulk_sectors.jl) |
| Address a coefficient by fusion trees and a junction | [Writable tree views](examples/03_tree_view.jl) · [Two bulk legs](examples/04_two_bulk_legs.jl) |
| Store boundary eigenvectors in a central tensor map | [Write and recover a right eigenstate](examples/05_state_columns.jl) |
