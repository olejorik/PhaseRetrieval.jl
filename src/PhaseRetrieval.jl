module PhaseRetrieval
using LinearAlgebra
using SampledDomains: CartesianDomain2D, dualRange
import SampledDomains: make_centered_domain2D

export SHdiversity
export psf, phwrap, subpsf
export Linear

include("types.jl")
include("utils.jl")
include("ShackHartmann/SHphase.jl")
include("PSF.jl")

"Wrap Phase" phwrap(x::Real) = mod(x, 2pi)
phwrap(NaN) = NaN



end
