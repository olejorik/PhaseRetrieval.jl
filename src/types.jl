

abstract type Domain{N}  end

"""
    PDplan(plan, diversity)

Construct `P`upil `D`iversity plan, which, if multiplied by the array `a` of proper dimensions, computes ``fft( diversity *a)``
"""
struct PDplan
    plan
    diversity::Array{<:Number}
end

PDplan(diversity::Array{<:Number}) = PDplan(plan_fft(diversity), diversity)

import Base: *
*(pd::PDplan, a::Array{<:Number}) = pd.plan * (pd.diversity .* a)
# export Base.:*

"""
`CameraChip(pixelsize = p, imagesize = (sizex, sizey), bitdepth = 8, channelbitdepth = 8` repesent a camera chip with a given pixel size, imagesize,
and bitdepth. If omitted,  bitdepth is set to 8.
"""
Base.@kwdef struct CameraChip
    pixelsize::Float64
    imagesize::Tuple{Int64, Int64}
    bitdepth::Int = 8
    channelbitdepth::Int = 8
end

#simple function to redefine the camera size
"""
    roi(cam::CameraChip, dims::Tuple)

Create camera that represents ROI of `cam`.
"""
# roi(cam::CameraChip, dims::Tuple) = CameraChip(cam.pixelsize, min.(cam.imagesize, dims),cam.bitdepth, cam.channelbitdepth) # correct approach if we want to limit roi to the hardware
roi(cam::CameraChip, dims::Tuple) = CameraChip(cam.pixelsize, dims,cam.bitdepth, cam.channelbitdepth) # sometimes this is needed
roi(cam::CameraChip, dims::Integer) = roi(cam, (dims, dims))

"""Fixed focal length lens"""
Base.@kwdef struct ImagingLens
    focallength::Float64
    aperture::Float64
end

"""
Imaging sensor or camera consists of a lens and camera chip.

Create a sensor using keyword syntaxis (in any order), `lens` and `cam` are required
```
ims = ImagingSensor(lens = ::ImagingLens, cam = ::CameraChip)
```

Other possible keywords to simulate misalingnment
```
    focal_distance (= lens.focallength  by default)
    lensorigin = [0.,0.] # aka nodal point, but expressed in length units
    # misalingment parameters (not realised yet, set to zero)
    α = 0.
    β = 0.
    γ = 0.
```
"""
Base.@kwdef struct ImagingSensor
    lens::ImagingLens
    cam::CameraChip
    focal_distance = lens.focallength # 
    lensorigin = [0.,0.] # aka nodal point, but expressed in length units
    # misalingment parameters (not realised yet, set to zero)
    α = 0.
    β = 0.
    γ = 0.
end

diaphragm(lens::ImagingLens, diam::Float64) = ImagingLens(lens.focallength, min.(lens.aperture, diam))
   
roi(ims::ImagingSensor, dims) = ImagingSensor(ims.lens, roi(ims.cam,dims),ims.focal_distance, ims.lensorigin,ims.α, ims.β, ims.γ)
# end


Base.@kwdef struct hwConfig
    cam::PhaseRetrieval.CameraChip
    f::Float64
    λ::Float64
    d::Float64
end

"""
    hwConfig(s::String, f, λ, d) creates a hardware configuration with
    `s` camera, lens with a focal length `f` and aperture `d` and using wavelenght
    `λ`.

## Example

    conf1 = hwConfig("UI1540", 300mm, 633nm,25mm) creates a configuration
    based on UI-1540 camera, with a 1 inch lens with focal length 300mm and He-Ne laser. 
"""
hwConfig(s::String, f, λ, d) = hwConfig(PhaseRetrieval.camerasdict[s], f, λ, d)

struct SimConfig
    name::String
    ims::PhaseRetrieval.ImagingSensor
    f::Float64
    λ::Float64
    d::Float64
    q::Int
    roi::CartesianDomain2D
    dualroi::CartesianDomain2D
    ap::Array{Float64,2}
    mask::Array{Float64,2}
    # diversity::Array{Float64,2} # not implemented
end

function SimConfig(name::String, ims::PhaseRetrieval.ImagingSensor, λ::Float64)
    q = PhaseRetrieval.upscaleFactor(ims, λ)
    f = ims.lens.focallength
    d = ims.lens.aperture
    roi = make_centered_domain2D(ims)
    dualroi = dualDomain(roi, q) * f * λ
    ap, mask = aperture(dualroi, d)
    return SimConfig(name, ims, f, λ, d, q, roi, dualroi, ap, mask)
end

function SimConfig(name::String, hw::hwConfig)
    ims = ImagingSensor(
        lens = ImagingLens(hw.f, hw.d),
        cam = hw.cam)
    λ = hw.λ
    return SimConfig(name, ims,  λ)
end


export hwConfig, SimConfig, ImagingSensor

