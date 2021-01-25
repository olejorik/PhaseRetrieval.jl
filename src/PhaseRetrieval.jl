module PhaseRetrieval
using LinearAlgebra
using FFTW
import Pkg
Pkg.add(url="https://github.com/olejorik/SampledDomains.jl")
# using SampledDomains
using SampledDomains: CartesianDomain2D, dualRange
import SampledDomains: make_centered_domain2D

export SHdiversity
export psf, phwrap, subpsf
export Linear

include("types.jl")
include("utils.jl")
include("ShackHartmann/SHphase.jl")
include("PSF.jl")





end
