
using FFTW
using MappedArrays

export AutoExposure, PSFMethods, PSFExposure
export wavelength, airysize, diversed_psfs, throughfocus, doflength




abstract type PSFExposure end
"""
    AutoExposure(s =1) <: PSFExposure

Type used to set exposure of the camera to `s`*m, where m is the exposure that images the maximal intensity of the frame as `1`. For instance, if the same scene is imaged with `AutoExposure(0.5)`, the brightest pixel will have value 0.5.
Values of `s` larger than 1 produce over-exposured images.

"""
struct AutoExposure <: PSFExposure
    scale::Float64
end


AutoExposure() = AutoExposure(1)





"""
    (cam::CameraChip)(
    intensity; exposure::PSFExposure=AutoExposure(), quantize=true, noise=(0, 0)
)

Calculate image generated by the camera for the given filed inensity.
    If used with complex array, first converts it to intesity
"""
function (cam::CameraChip)(
    imf::Array{T}; exposure::PSFExposure=AutoExposure(), quantize=true, noise=(0, 0)
) where {T<:Real}
    # imf = intensity(field)
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

(cam::CameraChip)(field::Array{Complex{T}}; kwargs...) where {T<:Real} =
    cam(collect(intensity(field)); kwargs...)

"""
    SimConfig(name::String, ims::PhaseRetrieval.ImagingSensor, λ::Float64; method = PSFMethods.Fourier()) creates forward-model
    simulation environment, using algorithm defined by `method`,
    for the PR problem obtained with  `ims` image sensor (camera + lens) at
    wavelength `λ`. `name` is a sting identifacator used for, for instace, plotting labels.

`SimConfig{PSFMethods.Fourier}` contains the following fields:

    - name::String
    - ims::PhaseRetrieval.ImagingSensor -- imaging sensor used in the configuration
    - λ::Float64 -- light wavelength
    - roi::CartesianDomain2D -- coordiantes in the image plane
    - dualroi::CartesianDomain2D -- coordinates in the pupil plane
    - ap::Array{Float64,2} -- amplitude of the pupil function
    - mask::Array{Float64,2} -- array of NaNs and 1s defining the aperture support
    - phases::Dict{String, PhaseBases.Phase} -- phase aberrrations present in the config
    - diversity::Dict{String, PhaseBases.Phase} -- phase diversities used for PSF generation
    - modulation::Dict{String, Array{Float64,2}} -- amplitude modulations used for calculation of the vectorial PSF


"""
struct SimConfig{PSFT<:PSFMethod}
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
    modulation::Dict{String,Array{Float64}}
    psfmethod::PSFT
end

function SimConfig(
    args...; method::T=default_forward_method, kwargs...
) where {T<:PSFMethod}
    return SimConfig{T}(args..., kwargs...)
end

function SimConfig{PSFMethods.Fourier}(name::String, ims::ImagingSensor, λ::Float64)
    q = upscaleFactor(ims, λ)
    f = ims.lens.focallength
    d = ims.lens.aperture
    roi = make_centered_domain2D(ims)
    dualroi = dualDomain(roi, q) * f * λ
    ap, mask = aperture(dualroi, d)
    return SimConfig{PSFMethods.Fourier}(
        name,
        ims,
        λ,
        roi,
        dualroi,
        ap,
        mask,
        Dict{String,Phase}(),
        Dict{String,Phase}(),
        Dict{String,Array{Float64,2}}(),
        PSFMethods.Fourier(),
    )
end

# Pretty printing
function show(io::IO, x::SimConfig)
    return print(
        io,
        """
        Simulation configuration with the properties:

        name:\t\t $(x.name)
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

algtype(x::SimConfig{T}) where {T} = T()

aperture(c::SimConfig) = error("Implement `aperture` method for $(typeof(c))")
aperture(c::SimConfig{PSFMethods.Fourier}) = c.ap

totalphase(c::SimConfig) = error("Implement `totalphase` method for $(typeof(c))")
totalphase(c::SimConfig{PSFMethods.Fourier}) =
    sum(collect, values(c.phases); init=zero(aperture(c)))

function pupilfield(c::SimConfig)
    return reduce(field, collect.(values(c.phases)); init=aperture(c))
end

function psf(c::SimConfig{T}; kwargs...) where {T}
    focalfield = toimageplane(pupilfield(c), algtype(c))
    return ret = c.ims.cam(focalfield; kwargs...)
    # TODO add noise
end

function diversed_psfs(c::SimConfig{T}; kwargs...) where {T}
    div_fields = vcat(
        [pupilfield(c)], [field(pupilfield(c), collect(d)) for d in values(c.diversity)]
    )
    return [c.ims.cam(toimageplane(f, algtype(c)); kwargs...) for f in div_fields]
end

focallength(c::SimConfig) = focallength(c.ims)
focaldistance(c::SimConfig) = focaldistance(c.ims)
wavelength(c::SimConfig) = c.λ
apdiameter(c::SimConfig) = apdiameter(c.ims)

(s::SimConfig)(phase) = psf(s.ap, collect(phase))

airysize(c::SimConfig) = 1.22 * wavelength(c) * focallength(c) / apdiameter(c)

import PhaseBases: ZernikeBW, ModalPhase
ZernikeBW(c::SimConfig{PSFMethods.Fourier}, order=10) =
    ZernikeBW(c.dualroi, apdiameter(c), order)

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
