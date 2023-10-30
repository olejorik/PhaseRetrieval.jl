module PhaseRetrieval
using LinearAlgebra
using Statistics
using FFTW
using SampledDomains: CartesianDomain2D, dualRange, dualDomain
import SampledDomains: make_centered_domain2D
import Base: show

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
export psf, subpsf, logrescale, airysize
export removepiston, removetiptilt, twinphase

export psf, SimConfig, AutoExposure, PSFmethod, Fourier, PSFExposure
export wavelength, airysize, diversed_psfs, throughfocus, doflength, focallength
export get_polarization_magnitudes, incoherent_psf
export FourierVectorial

# export Linear

include("types.jl")
include("utils.jl")
include("hardware.jl")

using .Hardware
export hwConfig, SimConfig, ImagingSensor, ImagingLens, CameraChip, roi, diaphragm
export camerasdict, lensesdict, m, mm, um, μm, nm
export focaldistance, focallength, apdiameter

include("ShackHartmann/SHphase.jl")

include("PRproblem.jl")

include("PSF.jl")
using .ForwardModel

export AbstractPRproblem, PRproblem, appsftoPR, solve, GS, GSparam, SimConfig

end
