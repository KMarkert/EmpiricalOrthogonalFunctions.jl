__precompile__()

module EmpiricalOrthoFuncs

using Statistics
using LinearAlgebra
using Logging
import Base.truncate

include("eoftypes.jl")
include("eoffuncs.jl")

export EmpiricalOrthoFunc,
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
