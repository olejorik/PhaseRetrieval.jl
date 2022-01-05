module PhaseRetrieval
# using LinearAlgebra
using Statistics
using FFTW
using SampledDomains: CartesianDomain2D, dualRange
import SampledDomains: make_centered_domain2D

import AlternatingProjections: Problem, Algorithm, solve, TwoSetsFP, IterativeAlgorithm, initial, tolerance, maxit, keephistory, snapshots # these we will need to change
using AlternatingProjections

export SHdiversity
export psf, phwrap, subpsf
export removepiston, twinphase

# export Linear

include("types.jl")
include("hardware.jl")
include("utils.jl")
include("ShackHartmann/SHphase.jl")
include("PSF.jl")

include("PRproblem.jl")

export AbstractPRproblem, PRproblem, appsftoPR, solve, GS, GSparam



end
