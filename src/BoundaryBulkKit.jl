module BoundaryBulkKit

import TensorKit
using TensorKit: Sector, Vect, block, dim, sectors, TensorSpace, ProductSpace, FusionTree,
    blockdim, blocksectors, codomain, dims, domain, fusiontrees, sectortype, SectorMismatch
using Dictionaries: Dictionary
using Strided: sreshape

export DrinfeldCenter, base_sector
export forget, halfbraiding_symbol, forget_vertex

export BoundaryBulkJunction, junctions
export CentralTensorMap, CentralTensorMapSpace

include("center.jl")
include("junction.jl")
include("central/space.jl")
include("central/tensor.jl")
include("central/views.jl")

end
