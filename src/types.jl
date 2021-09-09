
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


Base.@kwdef struct CameraChip
    pixelsize::Float64
    imagesize::Tuple{Int64, Int64}
end

#simple function to redefine the camera size
"""
    roi(cam::CameraChip, dims::Tuple)

Create camera that represents ROI of `cam`.
"""
roi(cam::CameraChip, dims::Tuple) = CameraChip(cam.pixelsize, min.(cam.imagesize, dims))
roi(cam::CameraChip, dims::Integer) = roi(cam, (dims, dims))

"""Fixed focal length lens"""
Base.@kwdef struct ImagingLens
    focallength::Float64
    aperture::Float64
end

"""
Imaging sensor or camera consists of a lens and camera chip.
"""
Base.@kwdef struct ImagingSensor
    lens::ImagingLens
    cam::CameraChip
    focal_distance = lens.focallength # 
    lensorigin = [0.,0.] # aka nodal point, but expressed in length units
    # misalingment parameters (not realised yet, set ot zero)
    α = 0.
    β = 0.
    γ = 0.
end

diaphragm(lens::ImagingLens, diam::Float64) = ImagingLens(lens.focallength, min.(lens.aperture, diam))
    
# end