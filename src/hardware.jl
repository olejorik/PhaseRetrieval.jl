# dictionaries with typical harware for ease of use with experimental data

const m = 1
const mm = 1e-3
const um = 1e-6
const μm = 1e-6
const nm = 1e-9


"""
    CameraChip(pixelsize, imagesize; <keyword args>)
    
Create a camera chip with a given square pixel size and image size.

You can also specify the bit depth of the camera and the channel bit
depth. For instance, 12-bit camera (1bitdepth =12`) often transfer its measurement packed
in 16-bit number (`channelbitdepth = 16`).

# Arguments
 - `pixelsize::Float64`: size of the pixel. Only square pixels are supported.
 - `imagesize::Tuple{Int, Int}`: (width, height) of the image.
 - `bitdepth::Int = 8`: bit detpth of the camera.
 - `channelbitdepth::Int = 8`: channel bit depth of the camera. 

Arguments can be specified in any order.

# Example
```julia
cam = CameraChip(;
    pixelsize=5.3um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8
)
```

"""
Base.@kwdef struct CameraChip
    pixelsize::Float64
    imagesize::Tuple{Int, Int}
    bitdepth::Int = 8
    channelbitdepth::Int = 8
end

#simple function to redefine the camera size
"""
    roi(cam, dims)

Create camera that represents centered regeion of interest (ROI) of `cam`.

ROI can be larger than the camera size! :-)

# Arguments

 - `cam`: camera (`CameraChip`) or imaging sensor (`ImagingSensor`)
 - `dims::Tuple{Int, Int}`: new (width, height) of the image. 
    Single number defines a suare ROI.
"""
roi(cam::CameraChip, dims::Tuple{Int, Int}) = CameraChip(cam.pixelsize, dims,cam.bitdepth, cam.channelbitdepth) # sometimes this is needed
roi(cam::CameraChip, dims::Integer) = roi(cam, (dims, dims))
# roi(cam::CameraChip, dims::Tuple) = CameraChip(cam.pixelsize, min.(cam.imagesize, dims),cam.bitdepth, cam.channelbitdepth) # correct approach if we want to limit roi to the hardware

"""
    ImagingLens(;focallength, aperture)

Create fixed focal length lens with focus = `focallength` and aperture dimater `aperture`.

Arguments can be speciefied in any order.

See also [`diaphragm`](@ref), [`lensesdict`](@ref).

# Example
```jldoctest
julia> lens = ImagingLens(300mm, 25mm)
ImagingLens(0.3, 0.025)
```
"""
Base.@kwdef struct ImagingLens
    focallength::Float64
    aperture::Float64
end

"""
    diaphragm(lens, diam)

Change the lens diameter of an `ImagingLens` or an `ImagingSensor`.

Diaphragm can be larger than the lens diameter (there is no "hardware limitation").

# Example
```jldoctest
julia> diaphragm(ImagingLens(300mm, 25mm), 10mm)
ImagingLens(0.3, 0.01)
```

"""
diaphragm(lens::ImagingLens, diam::Float64) = ImagingLens(lens.focallength, min.(lens.aperture, diam))


"""
    ImagingSensor(lens=lens, cam=cam)   

Imaging sensor or camera consisting of a lens `lens` and camera chip `cam`.

Create a sensor using keyword syntaxis (in any order), `lens` and `cam` are required
```
ims = ImagingSensor(lens = ::ImagingLens, cam = ::CameraChip)
```

See also [`roi`](@ref), [`diaphragm`](@ref).

# Arguments:
 - `lens::ImagingLens`: lens of the camera
 - `cam::CameraChip`: sensor chip of the camera

If lens and cam are strings, corresponding values from `lensesdict` and `camerasdict` are used.

Other possible keywords to simulate misalingnment (not implemented yet)

 - `focal_distance = lens.focallength`: Distance between the lens and the sensor.
 - `lensorigin = [0.,0.]`: nodal point, expressed in length units (`mm` or `um`)
 - Chip plane misalingment angles: 
    - `α = 0`: rotation in x
    - `β = 0`: rotation in y
    - `γ = 0`: rotation in z

# Example 

```jldoctest
julia> ImagingSensor(lens = ImagingLens(300mm, 25mm),
           cam = CameraChip(
               pixelsize=5.2um,
               imagesize=(1280, 1024),
               bitdepth=8,
               channelbitdepth=8
               )
           )
ImagingSensor(ImagingLens(0.3, 0.025), CameraChip(5.2e-6, (1280, 1024), 8, 8), 0.3, (0.0, 0.0), 0.0, 0.0, 0.0)
 
julia> diaphragm(ImagingSensor(lens = lensesdict["F300A25"], cam= camerasdict["UI1540"]), 10mm)
ImagingSensor(ImagingLens(0.3, 0.01), CameraChip(5.2e-6, (1280, 1024), 8, 8), 0.3, (0.0, 0.0), 0.0, 0.0, 0.0)
 
julia> ImagingSensor("F300A25",  "UI1540")
ImagingSensor(ImagingLens(0.3, 0.025), CameraChip(5.2e-6, (1280, 1024), 8, 8), 0.3, (0.0, 0.0), 0.0, 0.0, 0.0)
```
"""
Base.@kwdef struct ImagingSensor
    lens::ImagingLens
    cam::CameraChip
    focal_distance::Float64 = lens.focallength # 
    lensorigin::Tuple{Float64,Float64} = (0.,0.) # aka nodal point, but expressed in length units
    # misalingment parameters (not realised yet, set to zero)
    α::Float64 = 0.
    β::Float64 = 0.
    γ::Float64 = 0.
end


ImagingSensor(lensname::String, camname::String)=
    ImagingSensor(lens= lensesdict[lensname],  cam= camerasdict[camname])
# ImagingSensor(lensname::String, camname::String; kwargs...)=
#     ImagingSensor(lens= lensesdict[lensname],  cam= camerasdict[camname], kwargs...)



diaphragm(ims::ImagingSensor, diam) = ImagingSensor(diaphragm(ims.lens, diam), ims.cam,ims.focal_distance, ims.lensorigin,ims.α, ims.β, ims.γ)
   
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

"""
    SimConfig(name::String, ims::PhaseRetrieval.ImagingSensor, λ::Float64) creates forward-model 
    simulation environment for the PR problem obtained with  `ims` image sensor (camera + lens) at 
        wavelength `λ`. `name` is a sting identifacator used for, for instace, plotting labels.

`SimConfig` contains the following fields:
    ```julia
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
    ````

"""
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

(s::SimConfig)(phase) = psf(s.ap, phase)



"""
    `camerasdict` is a dictionary with often used cameras.

Use `keys(camerasdict)` to get the list of implemented cameras.
"""
camerasdict = Dict()

"""
    `lensesdict` is a dictionary with often used cameras.

Use `keys(lensesdict)` to get the list of implemented cameras.
"""
lensesdict = Dict()


# Cameras
camerasdict["UI1490"] = CameraChip(pixelsize = 1.67um, imagesize = (3840, 2748), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI1540"] = CameraChip(pixelsize = 5.2um, imagesize = (1280, 1024), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI1240"] = CameraChip(pixelsize = 5.3um, imagesize = (1280, 1024), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI2210"] = CameraChip(pixelsize = 9.9um, imagesize = (640, 480), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI3260"] = CameraChip(pixelsize = 5.86um, imagesize = (1936, 1216), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16
camerasdict["UI3860"] = CameraChip(pixelsize = 2.9um, imagesize = (1936, 1096), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16
camerasdict["MC203MG"] =CameraChip(pixelsize = 2.74um, imagesize = (4504, 4504), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16

# Lenses
lensesdict["F300A25"] = ImagingLens(300mm, 25mm)
lensesdict["F700A75"] = ImagingLens(700mm, 75mm)


export hwConfig, SimConfig, ImagingSensor, ImagingLens, CameraChip, roi, diaphragm
export camerasdict, lensesdict, m, mm, um, μm, nm