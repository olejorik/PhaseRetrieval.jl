module PhaseRetrieval

using LinearAlgebra
using Statistics
using FFTW
using SampledDomains: CartesianDomain2D, dualRange, dualDomain
using PhaseBases
using ImageCore
using AlternatingProjections

import Base: show, *
import SampledDomains: make_centered_domain2D, getranges

import AlternatingProjections: # these we will need to change
    Algorithm,
    initial,
    IterativeAlgorithm,
    keephistory,
    maxit,
    Problem,
    snapshots,
    solve,
    tolerance,
    TwoSetsFP

export
    # types
    AbstractPRproblem,
    AutoExposure,
    Fourier,
    FourierVectorial,
    GS,
    GSparam,
    PRproblem,
    PSFExposure,
    PSFmethod,
    SHdiversity,
    SimConfig,

    # methods
    airysize,
    appsftoPR,
    diversed_psfs,
    doflength,
    focallength,
    get_polarization_magnitudes,
    incoherent_psf,
    logrescale,
    psf,
    phwrap,
    removepiston,
    removetiptilt,
    solve,
    throughfocus,
    twinphase,
    wavelength

### Source files

# type system
include("types.jl")

# implementation helpers
include("utils.jl")

# generic functions
include("arrayutils.jl")
include("phasefunctions.jl")
include("psffunctions.jl")

#TODO refactor the rest
include("hardware.jl")

include("ShackHartmann/SHphase.jl")

# forward problem
include("PSF.jl")
include("VectorialPSF.jl")

# Inverse problem
include("PRproblem.jl")

end
