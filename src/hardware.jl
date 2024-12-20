export hwConfig, SimConfig, ImagingSensor, ImagingLens, CameraChip, roi, diaphragm, aperture
export camerasdict, lensesdict, m, mm, um, μm, nm
export focaldistance, focallength, apdiameter

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
depth. For instance, 12-bit camera (`bitdepth =12`) often transfer its measurement packed
in 16-bit number (`channelbitdepth = 16`).

If called on a complex field, output the image seen by the chip.

See also `camerasdict`, `storagetype`, `intensity`.

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

```jldoctest
julia> ap,_ = PhaseRetrieval.aperture(-1:.25:1,-1:.25:1,.8); ap
9×9 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> imf = PhaseRetrieval.toimageplane(ap);

julia> cam = camerasdict["UI3860"];

julia> cam(imf; exposure = AutoExposure(.7)) |> collect
9×9 Array{N4f12,2} with eltype FixedPointNumbers.N4f12:
 0.0051  0.0  0.0122  0.0427  0.0601  0.0427  0.0122  0.0  0.0051
 0.0     0.0  0.0     0.0     0.0     0.0     0.0     0.0  0.0
 0.0122  0.0  0.0286  0.1006  0.1411  0.1006  0.0286  0.0  0.0122
 0.0427  0.0  0.1006  0.3553  0.4987  0.3553  0.1006  0.0  0.0427
 0.0601  0.0  0.1411  0.4987  0.6999  0.4987  0.1411  0.0  0.0601
 0.0427  0.0  0.1006  0.3553  0.4987  0.3553  0.1006  0.0  0.0427
 0.0122  0.0  0.0286  0.1006  0.1411  0.1006  0.0286  0.0  0.0122
 0.0     0.0  0.0     0.0     0.0     0.0     0.0     0.0  0.0
 0.0051  0.0  0.0122  0.0427  0.0601  0.0427  0.0122  0.0  0.0051

```

"""
Base.@kwdef struct CameraChip
    pixelsize::Tuple{Float64,Float64}
    imagesize::Tuple{Int,Int}
    bitdepth::Int = 8
    channelbitdepth::Int = 8
end

CameraChip(pixelsize::Float64, imagesize, bitdepth, channelbitdepth) =
    CameraChip((pixelsize, pixelsize), imagesize, bitdepth, channelbitdepth)


#simple function to redefine the camera size
"""
    roi(cam, dims)

Create camera that represents centered regeion of interest (ROI) of `cam`.

ROI can be larger than the camera size! :-)

# Arguments

 - `cam`: camera (`CameraChip`) or imaging sensor (`ImagingSensor`)
 - `dims::Tuple{Int, Int}`: new (width, height) of the image.
    Single number defines a square ROI.
"""
function roi(cam::CameraChip, dims::Tuple{Int,Int})
    return CameraChip(cam.pixelsize, dims, cam.bitdepth, cam.channelbitdepth)
end # sometimes this is needed
roi(cam::CameraChip, dims::Integer) = roi(cam, (dims, dims))
# roi(cam::CameraChip, dims::Tuple) = CameraChip(cam.pixelsize, min.(cam.imagesize, dims),cam.bitdepth, cam.channelbitdepth) # correct approach if we want to limit roi to the hardware

"""
    storagetype[(bitdepth,channelbitdepth)]

Get the storage type returned by the combination of the bit depth and storage depth.
"""
const storagetype = (;
    :b8c8 => N0f8,
    :b8c16 => N8f8,
    :b10c16 => N6f10,
    :b12c16 => N4f12,
    :b14c16 => N2f14,
    :b16c16 => N0f16,
)

storagetypefun(::Val{8}, ::Val{8}) = N0f8
storagetypefun(::Val{8}, ::Val{16}) = N0f16
storagetypefun(::Val{10}, ::Val{16}) = N6f10
storagetypefun(::Val{12}, ::Val{16}) = N4f12
storagetypefun(::Val{14}, ::Val{16}) = N2f114
storagetypefun(::Val{16}, ::Val{16}) = N0f16

storagetypefun(::Val{-1}, ::Val{-1}) = Float64 # for abstract cameras without quantisation

"""
    getstoragetype(cam::CameraChip)

Get the storage type returned by the camera `cam`.
"""
function getstoragetype(cam::CameraChip)
    return storagetypefun(Val(cam.bitdepth), Val(cam.channelbitdepth))
end

"""
    ImagingLens(;focallength, aperture)

Create fixed focal length lens with focus = `focallength` and aperture dimater `aperture`.

Arguments can be speciefied in any order.

See also [`diaphragm`](@ref), [`lensesdict`](@ref).

# Example
```jldoctest
julia> lens = ImagingLens(300mm, 25mm)
Imaging Lens with
  f:	300.0 mm and
  D:	25.0 mm
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
Imaging Lens with
  f:	300.0 mm and
  D:	10.0 mm
```

"""
function diaphragm(lens::ImagingLens, diam::Float64)
    return ImagingLens(lens.focallength, min.(lens.aperture, diam))
end

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
Imaging Sensor made with:
  Imaging Lens with
  f:	300.0 mm and
  D:	25.0 mm and
  Camera with
	pixel size 5.2 × 5.2 μm²,
	(1280, 1024) frame and
	8/8 bit/channel output.

julia> diaphragm(ImagingSensor(lens = lensesdict["F300A25"], cam= camerasdict["UI1540"]), 10mm)
Imaging Sensor made with:
  Imaging Lens with
  f:	300.0 mm and
  D:	10.0 mm and
  Camera with
	pixel size 5.2 × 5.2 μm²,
	(1280, 1024) frame and
	8/8 bit/channel output.

julia> ImagingSensor("F300A25",  "UI1540")
Imaging Sensor made with:
  Imaging Lens with
  f:	300.0 mm and
  D:	25.0 mm and
  Camera with
	pixel size 5.2 × 5.2 μm²,
	(1280, 1024) frame and
	8/8 bit/channel output.
```

See also `focallength`, `focaldistance`, `apdiameter`.
"""
Base.@kwdef struct ImagingSensor
    lens::ImagingLens
    cam::CameraChip
    focal_distance::Float64 = lens.focallength #
    lensorigin::Tuple{Float64,Float64} = (0.0, 0.0) # aka nodal point, but expressed in length units
    # misalingment parameters (not realised yet, set to zero)
    α::Float64 = 0.0
    β::Float64 = 0.0
    γ::Float64 = 0.0
end

function ImagingSensor(lensname::String, camname::String)
    return ImagingSensor(; lens=lensesdict[lensname], cam=camerasdict[camname])
end
# ImagingSensor(lensname::String, camname::String; kwargs...)=
#     ImagingSensor(lens= lensesdict[lensname],  cam= camerasdict[camname], kwargs...)

# Setters
function diaphragm(ims::ImagingSensor, diam)
    return ImagingSensor(
        diaphragm(ims.lens, diam),
        ims.cam,
        ims.focal_distance,
        ims.lensorigin,
        ims.α,
        ims.β,
        ims.γ,
    )
end

function roi(ims::ImagingSensor, dims)
    return ImagingSensor(
        ims.lens,
        roi(ims.cam, dims),
        ims.focal_distance,
        ims.lensorigin,
        ims.α,
        ims.β,
        ims.γ,
    )
end
function focaldistance(ims::ImagingSensor, z::Float64)
    return ImagingSensor(ims.lens, ims.cam, z, ims.lensorigin, ims.α, ims.β, ims.γ)
end

# Getters
"""
    focaldistance(ims [,z])

Get/set the distance from the lens to the imaging plane.
"""
focaldistance(ims::ImagingSensor) = ims.focal_distance

focallength(ims::ImagingSensor) = focallength(ims.lens)
focallength(len::ImagingLens) = len.focallength
apdiameter(ims::ImagingSensor) = apdiameter(ims.lens)
apdiameter(len::ImagingLens) = len.aperture
numericalaperture(len::ImagingLens) = apdiameter(len) / focallength(len) / 2
numericalaperture(ims::ImagingSensor) = numericalaperture(ims.lens)

# Pretty printing
function show(io::IO, x::ImagingLens)
    return print(
        io,
        "Imaging Lens with \n  f:\t$(x.focallength/mm) mm and \n  D:\t$(x.aperture/mm) mm",
    )
end
function show(io::IO, x::ImagingSensor)
    return print(io, "Imaging Sensor made with:\n  $(x.lens) and \n  $(x.cam)")
end
function show(io::IO, x::CameraChip)
    return print(
        io,
        "Camera with \n\tpixel size $(x.pixelsize[1]/um) × $(x.pixelsize[2]/um) μm², \n\t$(x.imagesize) frame and \n\t$(x.bitdepth)/$(x.channelbitdepth) bit/channel output.",
    )
end

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
camerasdict["UI1490"] = CameraChip(;
    pixelsize=1.67um, imagesize=(3840, 2748), bitdepth=8, channelbitdepth=8
)
camerasdict["UI1540"] = CameraChip(;
    pixelsize=5.2um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8
)
camerasdict["UI1240"] = CameraChip(;
    pixelsize=5.3um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8
)
camerasdict["UI2210"] = CameraChip(;
    pixelsize=9.9um, imagesize=(640, 480), bitdepth=8, channelbitdepth=8
)
camerasdict["UI3260"] = CameraChip(;
    pixelsize=5.86um, imagesize=(1936, 1216), bitdepth=12, channelbitdepth=16
) #Bit depth of the images is 16
camerasdict["UI3860"] = CameraChip(;
    pixelsize=2.9um, imagesize=(1936, 1096), bitdepth=12, channelbitdepth=16
) #Bit depth of the images is 16
camerasdict["MC203MG"] = CameraChip(;
    pixelsize=2.74um, imagesize=(4504, 4504), bitdepth=12, channelbitdepth=16
) #Bit depth of the images is 16

# Lenses
lensesdict["F300A25"] = ImagingLens(300mm, 25mm)
lensesdict["F700A75"] = ImagingLens(700mm, 75mm)
