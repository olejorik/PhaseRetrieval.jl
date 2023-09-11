module PhaseRetrieval
using LinearAlgebra
using Statistics
using FFTW
using SampledDomains: CartesianDomain2D, dualRange, dualDomain
import Base: show
import SampledDomains: make_centered_domain2D, getranges
using PhaseBases
using ImageCore

import AlternatingProjections:
    Problem,
    Algorithm,
    solve,
    TwoSetsFP,
    IterativeAlgorithm,
    initial,
    tolerance,
    maxit,
    keephistory,
    snapshots # these we will need to change
using AlternatingProjections

export SHdiversity
export psf, phwrap, subpsf, logrescale
export removepiston, removetiptilt, twinphase

# export Linear

include("types.jl")
include("hardware.jl")
include("utils.jl")
include("ShackHartmann/SHphase.jl")
include("PSF.jl")
include("VectorialPSF.jl")

include("PRproblem.jl")

export AbstractPRproblem, PRproblem, appsftoPR, solve, GS, GSparam

end
