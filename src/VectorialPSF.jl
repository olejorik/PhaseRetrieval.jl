function get_cosines(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    σx, σy = getranges(ddom)
    # σs = (σz(σ...) - 1 for σ in ddom)
    σs = (σz(σ...) for σ in ddom)
    return (σx=σx, σy=σy, σz=σs)
end

function get_polarization_magnitudes(ddom::CartesianDomain2D, inaperture::BitArray)
    # σx, σy = getranges(ddom)
    exx = Array{Float64}(undef, size(ddom))
    eyx = Array{Float64}(undef, size(ddom))
    ezx = Array{Float64}(undef, size(ddom))
    exy = Array{Float64}(undef, size(ddom))
    eyy = Array{Float64}(undef, size(ddom))
    ezy = Array{Float64}(undef, size(ddom))
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

function get_polarization_magnitudes(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    inaperture = (!isnan).(conf.mask)
    return get_polarization_magnitudes(ddom, inaperture)
end

struct FourierVectorial <: PSFmethod
    fplan::FFTW.FFTWPlan{ComplexF64,-1,false}
    modulations
end

function FourierVectorial(dom::CartesianDomain2D, mask::Array)
    pf = FFTW.plan_fft(Array{ComplexF64}(undef, size(dom)))
    inaperture = (!isnan).(mask)
    return FourierVectorial(pf, get_polarization_magnitudes(dom, inaperture))
end

function SimConfig{FourierVectorial}(name::String, ims::ImagingSensor, λ::Float64)
    q = upscaleFactor(ims, λ)
    f = ims.lens.focallength
    d = ims.lens.aperture
    roi = make_centered_domain2D(ims)
    dualroi = dualDomain(roi, q) * f * λ
    ap, mask = aperture(dualroi, d)
    alg = FourierVectorial(dualroi, mask)
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

function toimageplane(field, alg::FourierVectorial, pol_coeffs)
    return [
        fftshift(alg.fplan * ifftshift((p * m) .* field)) for
        (m, p) in zip(alg.modulations, pol_coeffs)
    ]
end

function toimageplane(field, alg::FourierVectorial)
    return toimageplane(field, alg::FourierVectorial, ones(6))
end

aperture(c::SimConfig{FourierVectorial}) = c.ap

export FourierVectorial
