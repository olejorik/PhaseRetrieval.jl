
struct PRForwardModel
    name::String
    ims::PhaseRetrieval.ImagingSensor
    f::Float64
    λ::Float64
    d::Float64
    q::Int
    roi::CartesianDomain2D
    dualroi::CartesianDomain2D
    ap::Array{Float64,2}
    phasepixel::Array{Float64,2}
    phaseZer::Array{Float64,2}
    mask::Array{Float64,2}
    diversity::Vector{Array{Float64,2}}
    noise::Float64
end

function PRForwardModel(
    name::String, ims::PhaseRetrieval.ImagingSensor; λ=633nm, diversity=[], noise=0
)
    q = PhaseRetrieval.upscaleFactor(ims, λ)
    if q > 1
        @info "Camera undersamples the PSF, using upscaling with factor $q"
    end
    f = ims.lens.focallength
    d = ims.lens.aperture
    roi = make_centered_domain2D(ims)
    dualroi = dualDomain(roi, q) * f * λ
    ap, mask = aperture(dualroi, d)
    phasepixel = zero(ap)
    phaseZer = zero(ap)
    return PRForwardModel(
        name,
        ims,
        f,
        λ,
        d,
        q,
        roi,
        dualroi,
        ap,
        phasepixel,
        phaseZer,
        mask,
        diversity,
        noise,
    )
end
