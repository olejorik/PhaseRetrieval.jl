# dictionaries with typical harware for ease of use with experimental data

const m = 1
const mm = 10e-3
const um = 10e-6
const μm = 10e-6
const nm = 10-9

camerasdict = Dict()
lensesdict = Dict()


# Cameras
camerasdict["UI1490"] = PhaseRetrieval.CameraChip(pixelsize = 1.67um, imagesize = (3840, 2748), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI1540"] = PhaseRetrieval.CameraChip(pixelsize = 5.2um, imagesize = (1280, 1024), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI2210"] = PhaseRetrieval.CameraChip(pixelsize = 9.9um, imagesize = (640, 480), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI3260"] = PhaseRetrieval.CameraChip(pixelsize = 5.86um, imagesize = (1936, 1216), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16
camerasdict["UI3860"] = PhaseRetrieval.CameraChip(pixelsize = 2.9um, imagesize = (1936, 1096), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16
camerasdict["MC203MG"] = PhaseRetrieval.CameraChip(pixelsize = 2.74um, imagesize = (4504, 4504), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16

# Lenses
lensesdict["F300A25"] = PhaseRetrieval.ImagingLens(300mm, 25mm)
lensesdict["F700A75"] = PhaseRetrieval.ImagingLens(700mm, 75mm)

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

export hwConfig, SimConfig, m, mm, um, μm, nm