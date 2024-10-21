export incoherent_psf, vectorial_psf, vectorial_diversed_psfs

"""
    σz(σx, σy)

Calculate σz from the values of the other cosines
"""
σz(σx, σy) = sqrt(1 + 0im - σx^2 - σy^2) # convert to complex to allow complex sigmas

function get_cosines(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    σx, σy = getranges(ddom)
    σs = (σz(σ...) for σ in ddom)
    return (σx=σx, σy=σy, σz=σs)
end

function get_polarization_magnitudes(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    inaperture = (!isnan).(conf.mask)
    # σx, σy = getranges(ddom)
    # exx = Array{Float64}(undef, size(ddom))
    # eyx = Array{Float64}(undef, size(ddom))
    # ezx = Array{Float64}(undef, size(ddom))
    # exy = Array{Float64}(undef, size(ddom))
    # eyy = Array{Float64}(undef, size(ddom))
    # ezy = Array{Float64}(undef, size(ddom))
    exx = zeros(size(ddom))
    eyx = zeros(size(ddom))
    ezx = zeros(size(ddom))
    exy = zeros(size(ddom))
    eyy = zeros(size(ddom))
    ezy = zeros(size(ddom))
    for (i, (y, x)) in enumerate(ddom)
        if inaperture[i]
            zfactor = 1 / (1 + σz(x, y))
            exx[i] = 1 - x^2 * zfactor
            eyx[i] = -x * y * zfactor
            ezx[i] = -x
            exy[i] = -x * y * zfactor
            eyy[i] = 1 - y^2 * zfactor
            ezy[i] = -y
        end
    end
    return (exx=exx, eyx=eyx, ezx=ezx, exy=exy, eyy=eyy, ezy=ezy)
end

function get_obliquity_factor(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    inaperture = (!isnan).(conf.mask)
    ret = ones(size(ddom))
    for (i, (y, x)) in enumerate(ddom)
        if inaperture[i]
            ret[i] = 1 / sqrt(real(σz(x, y)))
        end
    end
    return ret
end

"""
    incoherent_psf(array<:Array)

Calculate 2D array corresponding to the incoherent sum of the PSF formed by slices `array[:,:,ind...]` for all possible indici `ind`.

```jldoctest
julia> incoherent_psf(ones(3,3,2,3))
3×3 Matrix{Float64}:
 0.0    0.0  0.0
 0.0  486.0  0.0
 0.0    0.0  0.0

julia> incoherent_psf(ones(4,4,2,2,1))
4×4 Matrix{Float64}:
 0.0  0.0     0.0  0.0
 0.0  0.0     0.0  0.0
 0.0  0.0  1024.0  0.0
 0.0  0.0     0.0  0.0

julia> incoherent_psf(ones(4,4))
4×4 Matrix{Float64}:
 0.0  0.0    0.0  0.0
 0.0  0.0    0.0  0.0
 0.0  0.0  256.0  0.0
 0.0  0.0    0.0  0.0

```
"""
function incoherent_psf(array::Array)
    arr_dims = ndims(array)
    arr_dims > 1 || error("Array should have at least 2 dimensions")

    ret = zeros(size(array)[1:2])
    for f in eachslice(array; dims=Tuple(3:arr_dims))
        ret .+= psf(f)
    end

    return ret
end

function _vectorial_sim_config(
    name::String, ims::ImagingSensor, λ::Float64, ::RandomPolarization
)
    conf = SimConfig(name, ims, λ)
    o = get_obliquity_factor(conf)
    pol_modes = get_polarization_magnitudes(conf)
    # for k in keys(pol_modes)
    #     conf.modulation[String(k)] = pol_modes[k] .* o
    # end
    conf.modulation["vectorial"] = stack(values(pol_modes)) .* o
    return conf
end

function _vectorial_sim_config(
    name::String, ims::ImagingSensor, λ::Float64, p::LinearPolarization
)
    conf = SimConfig(name, ims, λ)
    pol_modes = get_polarization_magnitudes(conf)
    o = get_obliquity_factor(conf)
    a = p.orientation_angle
    # for c in ["x", "y", "z"]
    #     conf.modulation[c] =
    #         (cos(a) * pol_modes[Symbol("e$(c)x")] + sin(a) * pol_modes[Symbol("e$(c)y")]) .*
    #         o
    # end
    conf.modulation["vectorial"] =
        stack([
            (cos(a) * pol_modes[Symbol("e$(c)x")] + sin(a) * pol_modes[Symbol("e$(c)y")])
            for c in ["x", "y", "z"]
        ]) .* o
    return conf
end

function vectorial_psf(c::SimConfig{T}; kwargs...) where {T}
    # focalfield = toimageplane(pupilfield(c), algtype(c))
    # return ret = c.ims.cam(focalfield, exposure, quantize, noise)
    # TODO add noise
    # TODO make consistent with the scalar
    haskey(c.modulation, "vectorial") ||
        error("Config named `$(c.name)` is not initialised for vectorial simulation")
    return c.ims.cam(incoherent_psf(pupilfield(c) .* c.modulation["vectorial"]); kwargs...)
end

function vectorial_diversed_psfs(c::SimConfig{T}; kwargs...) where {T}
    haskey(c.modulation, "vectorial") ||
        error("Config named `$(c.name)` is not initialised for vectorial simulation")

    div_fields =
    # vcat(
    # [pupilfield(c)],
    [cis.(d) for d in values(c.diversity)]
    # )
    # return [
    #     c.ims.cam(toimageplane(f, algtype(c)), exposure, quantize, noise) for
    #     f in div_fields
    # ]
    return [
        c.ims.cam(
            incoherent_psf(pupilfield(c) .* f .* c.modulation["vectorial"]); kwargs...
        ) for f in div_fields
    ]
end
