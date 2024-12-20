module PhaseRetrieval
using LinearAlgebra
using Statistics
using FFTW
using SampledDomains: CartesianDomain2D, dualRange, dualDomain
import Base: show, collect
import SampledDomains: make_centered_domain2D, getranges
using PhaseBases
using ImageCore

import AlternatingProjections:
    Problem,
    Algorithm,
    solve,
    TwoSetsFP,
    IterativeAlgorithm,
    APparam,
    initial,
    tolerance,
    maxit,
    keephistory,
    snapshots # these we will need to change
using AlternatingProjections

export SHdiversity
export psf, subpsf, logrescale
export removepiston, removetiptilt, twinphase
export gaussian_apodization

# export Linear

include("types.jl")
include("Field.jl")
include("hardware.jl")
include("utils.jl")
include("functions.jl")
include("methods.jl")
include("ShackHartmann/SHphase.jl")
include("PSF.jl")
include("VectorialPSF.jl")

include("PRproblem.jl")

export AbstractPRproblem,
    PRproblem,
    PRproblemSat,
    PDPRproblem,
    prproblem,
    appsftoPR,
    solve,
    stochasticSolveKeep,
    stochasticSolve,
    GS,
    GSparam,
    AP,
    DRAP,
    APparam,
    DRAPparam

end
