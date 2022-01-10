__precompile__()

module EmpiricalOrthogonalFunctions

using Statistics
using LinearAlgebra
using Logging
import Base.truncate

include("eoftypes.jl")
include("eoffuncs.jl")

export EmpiricalOrthogonalFunction,
    pcs,
    eofs,
    correlationmap,
    covariancemap,
    eigenvalues,
    variancefraction,
    totalanomalyvar,
    northtest,
    projectfield,
    reconstruct,
    orthorotation,
    orthorotation!,
    truncate,
    truncate!

end # module
