using FFTW
using MappedArrays

export psf, SimConfig, AutoExposure, PSFmethod, Fourier, PSFExposure
export wavelength, airysize, diversed_psfs, throughfocus, doflength, focallength
export get_polarization_magnitudes, incoherent_psf, intensity, toimageplane

abstract type PSFmethod end

"""
    Fourier(A::Array)

FFT-based method of simulation of a PSF.

See also `psf`, `SimConfig`.
"""
struct Fourier <: PSFmethod
    fplan::FFTW.FFTWPlan{ComplexF64,-1,false}
end

Fourier(A::Array) = Fourier(FFTW.plan_fft(A))

struct ChirpZ <: PSFmethod end
struct MVM <: PSFmethod end

const default_forward_method = Fourier

abstract type PSFExposure end
struct AutoExposure <: PSFExposure
    scale::Float64
end

struct GainExposure <: PSFExposure
    scale::Float64
end

AutoExposure() = AutoExposure(1)

myabs2(x) = abs(x)^2
# myabs2(x) = abs2(x)

field(amplitude, phase) = amplitude .* exp.(1im * collect(phase)) # use collect here for ModalPhase

"""
    psf(amplitude, phase) -> psfimage
    psf(pupilfield) -> psfimage


Calculate psf for given amplitude and phase. Both PSF and pupil arrays are assumed
    to have the origin in their center pixel (`(N+1)÷2` for dimension size `N`).

Examples
====
```jldoctest
julia> ap,_ = PhaseRetrieval.aperture(-1:.2:1,-1:.2:1,.8); ap
11×11 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> Int.(round.(psf(ap)))
11×11 Matrix{Int64}:
 18  3   1   2    0    1    0   2   1  3  18
  3  0   6   5    0    1    0   5   6  0   3
  1  6  11   4    1    5    1   4  11  6   1
  2  5   4   1   20   38   20   1   4  5   2
  0  0   1  20   78  118   78  20   1  0   0
  1  1   5  38  118  169  118  38   5  1   1
  0  0   1  20   78  118   78  20   1  0   0
  2  5   4   1   20   38   20   1   4  5   2
  1  6  11   4    1    5    1   4  11  6   1
  3  0   6   5    0    1    0   5   6  0   3
 18  3   1   2    0    1    0   2   1  3  18

```

```jldoctest
julia> psf(ones(3,3),zeros(3,3))
3×3 Matrix{Float64}:
 0.0   0.0  0.0
 0.0  81.0  0.0
 0.0   0.0  0.0

julia> psf(ones(Complex{Float64},4,4))
4×4 Matrix{Float64}:
 0.0  0.0    0.0  0.0
 0.0  0.0    0.0  0.0
 0.0  0.0  256.0  0.0
 0.0  0.0    0.0  0.0

```
"""
psf(field::AbstractArray) = myabs2.(toimageplane(field))
# psf(field::AbstractArray) = abs.(toimageplane(field)) .^2

psf(amplitude, phase) = psf(field(amplitude, phase))

subpsf(field::AbstractArray, Q::Integer) = psf(subdivide_sum(field, Q))
subpsf(amplitude, phase, Q::Integer) = psf(subdivide_sum(field(amplitude, phase), Q))

"""
    toimageplane(field [, method = Fourier()])

Calculate the field in the focal plane from the field in the pupil plane using `method`.

Available methods:
 - `Fourier()`
 - `ChirpZ`  -- not implemented
 - `MVM` -- not implemented
"""
toimageplane(field, alg::Fourier) = fftshift(alg.fplan * ifftshift(field))
toimageplane(field) = toimageplane(field, Fourier(field))

"""
    intensity(field)

Calculate intensity of the field.

# Examples

```jldoctest
julia> ap,_ = PhaseRetrieval.aperture(-1:.2:1,-1:.2:1,.8);

julia> imf = PhaseRetrieval.toimageplane(ap);

julia> PhaseRetrieval.intensity(imf) .|> round
11×11 Matrix{Float64}:
 18.0  3.0   1.0   2.0    0.0    1.0    0.0   2.0   1.0  3.0  18.0
  3.0  0.0   6.0   5.0    0.0    1.0    0.0   5.0   6.0  0.0   3.0
  1.0  6.0  11.0   4.0    1.0    5.0    1.0   4.0  11.0  6.0   1.0
  2.0  5.0   4.0   1.0   20.0   38.0   20.0   1.0   4.0  5.0   2.0
  0.0  0.0   1.0  20.0   78.0  118.0   78.0  20.0   1.0  0.0   0.0
  1.0  1.0   5.0  38.0  118.0  169.0  118.0  38.0   5.0  1.0   1.0
  0.0  0.0   1.0  20.0   78.0  118.0   78.0  20.0   1.0  0.0   0.0
  2.0  5.0   4.0   1.0   20.0   38.0   20.0   1.0   4.0  5.0   2.0
  1.0  6.0  11.0   4.0    1.0    5.0    1.0   4.0  11.0  6.0   1.0
  3.0  0.0   6.0   5.0    0.0    1.0    0.0   5.0   6.0  0.0   3.0
 18.0  3.0   1.0   2.0    0.0    1.0    0.0   2.0   1.0  3.0  18.0


```

"""
intensity(field) = mappedarray(myabs2, field)
intensity(fields::Vector{T}) where {T<:Array} = sum(intensity, fields)

function (cam::CameraChip)(
    field, exposure::PSFExposure=AutoExposure(), quantize=true, noise=(0, 0)
)
    imf = intensity(field)
    maxint = maximum(imf) / exposure.scale
    if noise != (0, 0) # TODO add docs that noise is scaled to the smallest bit
        imf = mappedarray(
            x -> x + maxint / (2^cam.bitdepth) * (randn() * noise[1] + noise[2]), imf
        )
    end
    storagetype = quantize ? getstoragetype(cam) : Float64
    # storagetype = N4f12 #debug -- it's three times faster compared with the unnknown type
    return mappedarray(scaleminmax(storagetype, 0, maxint), imf)
end

"""
    SimConfig(name::String, ims::PhaseRetrieval.ImagingSensor, λ::Float64; method = Fourier()) creates forward-model
    simulation environment, using algorithm defined by `method`,
    for the PR problem obtained with  `ims` image sensor (camera + lens) at
    wavelength `λ`. `name` is a sting identifacator used for, for instace, plotting labels.

`SimConfig{Fourier}` contains the following fields:

    - name::String
    - ims::PhaseRetrieval.ImagingSensor -- imaging sensor used in the configuration
    - λ::Float64 -- light wavelength
    - roi::CartesianDomain2D -- coordiantes in the image plane
    - dualroi::CartesianDomain2D -- coordinates in the pupil plane
    - ap::Array{Float64,2} -- amplitude of the pupil function
    - mask::Array{Float64,2} -- array of NaNs and 1s defining the aperture support
    - phases::Dict{String, PhaseBases.Phase} -- phase aberrrations present in the config
    - diversity::Dict{String, PhaseBases.Phase} -- phase diversities used for PSF generation
    - alg::PSFmethod -- method used for PSF calculation


"""
struct SimConfig{T<:PSFmethod}
    name::String
    ims::ImagingSensor
    λ::Float64
    # q::Int
    roi::CartesianDomain2D
    dualroi::CartesianDomain2D
    ap::Array{Float64,2}
    mask::Array{Float64,2}
    phases::Dict{String,Phase}
    diversity::Dict{String,Phase} # TODO add defocus calculation from focaldistance
    alg::T
end

# function SimConfig(
#     args...; method::T=default_forward_method, kwargs...
# ) where {T<:PSFmethod}
#     return SimConfig{T}(args..., kwargs...)
# end

function SimConfig(name::String, ims::ImagingSensor, λ::Float64)
    q = upscaleFactor(ims, λ)
    f = ims.lens.focallength
    d = ims.lens.aperture
    roi = make_centered_domain2D(ims)
    dualroi = dualDomain(roi, q) * f * λ
    ap, mask = aperture(dualroi, d)
    alg = Fourier(ap)
    return SimConfig(
        name,
        ims,
        λ,
        roi,
        dualroi,
        ap,
        mask,
        Dict{String,Phase}(),
        Dict{String,Phase}(),
        alg,
    )
end

# Pretty printing
function show(io::IO, x::SimConfig)
    return print(
        io,
        """
        Simulation configuration with the properties:

        name:\t\t $(x.name)
        algorithm:\t $(typeof(x.alg))
        wavelength:\t $(x.λ)
        $(x.ims)
        with wavefront aberrrations: $(keys(x.phases))
        with phase diversities $(keys(x.diversity))
        """,
        # "PSF simulation configuration with the properties:

        # with wavefront aberrrations: $(x.phases)\n
        # with phase diversities ($x.diversity)\n",
    )
end

algtype(x::SimConfig{T}) where {T} = T

aperture(c::SimConfig) = error("Implement `aperture` method for $(typeof(c))")
aperture(c::SimConfig{Fourier}) = c.ap

totalphase(c::SimConfig) = error("Implement `totalphase` method for $(typeof(c))")
totalphase(c::SimConfig{Fourier}) = sum(collect, values(c.phases); init=zero(aperture(c)))

function pupilfield(c::SimConfig)
    return reduce(field, collect.(values(c.phases)); init=aperture(c))
end

function psf(
    c::SimConfig{T}; noise=(1, 0), exposure=AutoExposure(), quantize=true
) where {T}
    focalfield = toimageplane(pupilfield(c), c.alg)
    return ret = c.ims.cam(focalfield, exposure, quantize, noise)
    # TODO add noise
end

function diversed_psfs(
    c::SimConfig{T}; noise=(0, 0), exposure=AutoExposure(), quantize=true
) where {T}
    div_fields = vcat(
        [pupilfield(c)], [field(pupilfield(c), collect(d)) for d in values(c.diversity)]
    )
    return [
        c.ims.cam(toimageplane(f, c.alg), exposure, quantize, noise) for f in div_fields
    ]
end

"""
    incoherent_psf(
    c::SimConfig{T},
    amp_diversity::Vector;
    noise=(1, 0),
    exposure=AutoExposure(),
    quantize=true,
) where {T}

Calculate incoherent sum of PSFs formed by amplitude diversities.
"""
function incoherent_psf(
    c::SimConfig{T},
    amp_diversity::Vector;
    noise=(1, 0),
    exposure=AutoExposure(),
    quantize=true,
) where {T}
    ret = zero(aperture(c))
    focalfields = [toimageplane(a .* pupilfield(c), c.alg) for a in amp_diversity]

    ret = c.ims.cam(focalfields, exposure, quantize, noise)
    # TODO add noise

    return collect(c.ims.cam(focalfields, exposure, quantize, noise))
end

focallength(c::SimConfig) = focallength(c.ims)
focaldistance(c::SimConfig) = focaldistance(c.ims)
wavelength(c::SimConfig) = c.λ
apdiameter(c::SimConfig) = apdiameter(c.ims)

(s::SimConfig)(phase) = psf(s.ap, collect(phase))

airysize(c::SimConfig) = 1.22 * wavelength(c) * focallength(c) / apdiameter(c)

import PhaseBases: ZernikeBW, ModalPhase
ZernikeBW(c::SimConfig{Fourier}, order=10) = ZernikeBW(c.dualroi, apdiameter(c), order)

(ph::Phase)(c::SimConfig) = set_phase!(c, ph, "aberration")
(ph::Phase)(c::SimConfig, phasenature::String) = set_phase!(c, ph, phasenature)
set_phase!(c::SimConfig, ph::Phase, phasenature::String) = (c.phases[phasenature] = ph; c)

# TODO think how to restructure all this
# utils add-on
function ap_ratio(c::SimConfig)
    return ap_ratio(c.ims, c.λ)
end

function upscaleFactor(c::SimConfig)
    return upscaleFactor(c.ims, c.λ)
end

"""
    σz(σx, σy)

Calculate σz from the values of the other cosines
"""
σz(σx, σy) = sqrt(1 - σx^2 - σy^2)

"""
    throughfocus(conf::SimConfig, Δz)

Calculate the phase aberration (defocus) corresponding to the axial disaplcament on length
`Δz`.
"""
function throughfocus(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    return collect(σz(σ...) - 1 for σ in ddom)
end

throughfocus(conf::SimConfig, Δz) = throughfocus(conf::SimConfig) .* (Δz * 2π / conf.λ)
doflength(conf::SimConfig) = conf.λ / (numericalaperture(conf.ims)^2)
