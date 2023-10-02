# ```@meta
# CurrentModule = PhaseRetrieval
# DocTestSetup = quote
#     using PhaseRetrieval
# end
# ```

# # Hardware types

# ## Camera chip
# [`CameraChip`](@ref) represents a device that outputs sampled field intensity quantized
# and packed in 8- or 16-bit integer per color channel (currently, only monochrome
# cameras are implemented).

# This code sets up a camera with 5.2μm pixel, 1280 × 1024 frame, `bitdepth` and `channeldepth`
#  of 8 bits:

using PhaseRetrieval
using PhasePlots
using CairoMakie
using ImageCore
cam = CameraChip(; pixelsize=5.2um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8)

# Cameras can be also created using [`camerasdict`](@ref) dictionary.
# Let's consider examples of 8-bit and 12-bit cameras:
cam8bit = camerasdict["UI1240"]

#  and
cam12bit = camerasdict["UI3860"]

#  and compare the results of using these two cameras to measure this
#  field created by a circular aperture:
ap, _ = PhaseRetrieval.aperture(-1:0.025:1, -1:0.025:1, 0.9)

showarray(ap)

# The field itself is represented by a complex array, which has this amplitude and phase (we
# can see that due to the aperture symmetry, the field in the focal plane is real, actually,
# and the phase take only values 0 and π):
imf = PhaseRetrieval.toimageplane(ap)
showarray(abs.(imf))
#
showphasetight(angle.(imf))[1]

# And now we can quantize this field by our 8-bit and 12-bit example cameras
collect(cam8bit(imf))
#
collect(cam12bit(imf))

# !!! note
#     We have used `collect` here, as `CameraChip` returns a mapped view of the field.

# We can see that for our 12-bit camera, the measurements are represented by the last
# 12 bits of a 16-bit number
# ([`N4f12`](https://github.com/JuliaMath/FixedPointNumbers.jl)).

# We can hardly see the difference in the camera outputs. Compare 8bit (left)
# with 12 bit (right):
Gray.(mosaicview(cam8bit(imf), cam12bit(imf); nrow=1, npad=5, fillvalue=1))

# The difference is visible in the logarithmic scale. compare 8bit:
#  and 12 bit
Gray.(
    mosaicview(
        logrescale(float.(cam8bit(imf))),
        logrescale(float.(cam12bit(imf)));
        nrow=1,
        npad=5,
        fillvalue=1,
    )
)

# The difference appears because of the thresholding done by the quantization.
# Without the quantizations, the fields are the same:
Gray.(
    mosaicview(
        logrescale(float.(cam8bit(imf, AutoExposure(1), false))),
        logrescale(float.(cam12bit(imf, AutoExposure(1), false)));
        nrow=1,
        npad=5,
        fillvalue=1,
    )
)

# We can, however, boost exposure 16 (=2⁴) times of the 8bit camera to see more rings:
Gray.(PhaseRetrieval.logrescale(float.(cam8bit(imf, AutoExposure(16)))))
# We lost, of course, the information in the center of the PSF due to
# the oversaturation of the camera pixels.

# !!! note
#     Here we have used already sampled optical filed `imf`, so there is no difference in
#     the sampling rate performed by the cameras. This effect will be demonstrated in the next example.

# ## SimConfig
# `SimConfig` contains the necessary information for simulations.

# Let's set up a simulation environment matching the following hardware setup:
# a beam with a footprint of 1 inch (25 mm) diameter is focused with a lens of 300 mm
# focal length, and the PSF is registered with [UI-1240 camera](https://en.ids-imaging.com/store/products/cameras/ui-1240le.html).
# For both lens and camera, we can use structures with self-explanatory names
# [`ImagingLens`](@ref) and [`CameraChip`](@ref), which we combine in one structure called [`ImagingSensor`](@ref):

using PhaseRetrieval
lens = ImagingLens(300mm, 25mm)
cam = CameraChip(; pixelsize=5.3um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8)
ims = ImagingSensor(; lens=lens, cam=cam)

# Now we can save all these definitions in a simulation config [`SimConfig`](@ref). We also specify the wavelength here:

conf1 = SimConfig("full_aperture", ims, 633nm)

# This creates an aperture array of correct dimensions which is suitable for the generation of a PSF using Fourier methods.
# If we check near the central pixel, we'll see that for this configuration the PSF is almost one pixel wide

p = psf(conf1.ap)
using CairoMakie # hide
CairoMakie.activate!(; type="png") # hide
heatmap(rotr90(p[503:523, 631:651]); axis=(aspect=DataAspect(),))

# Indeed, the Airy pattern should be about 9 microns wide

print("Airy size is 1.22λ/NA = $(airysize(conf1)/ um) μm")

# We might thus want to consider a smaller numerical aperture:

lens2 = PhaseRetrieval.diaphragm(lens, 10mm)
ims2 = PhaseRetrieval.ImagingSensor(; lens=lens2, cam=cam)
conf2 = SimConfig("10mm aperture", ims2, 633nm)
p2 = psf(conf2.ap)
heatmap(rotr90(p2[503:523, 631:651]); axis=(aspect=DataAspect(),))

# #### Faster creation of an `ImagingSensor`
# Some often used cameras are saved in [`camerasdict`](@ref) and [`lensesdict`](@ref) dictionaries

keys(camerasdict)

# So the imaging sensor can be created as

ims = ImagingSensor(; cam=camerasdict["UI1240"], lens=lensesdict["F300A25"])

#
# `SimConfig` contains the necessary information for simulations.

fieldnames(typeof(conf2))

# For instance, it contains the aperture mask.

using PhasePlots
showarray(conf2.ap)

# The dimensions of the mask correspond to the dimensions of the sampled image plane,
# but the overall size corresponds to the inverse of the pixel size.
# This information is contained in `dualroi` field and can be used to construct the Zernike basis.

conf2.dualroi

# This can be used to create the Zernike basis

using PhaseBases
basis = ZernikeBW(conf2.dualroi, apdiameter(conf2), 10);
showphase(basis.elements[15] .* conf2.mask)[1]

# Or the same picture without unnecessary information (by default all phases will be shown scaled to (-π. π])

showphasetight(basis.elements[15] .* conf2.mask)[1]

# This is a combination of some low-order Zernike polynomials

phase = ModalPhase([4, 6, 15, 16], [2, 1, 0.4, 0.3] * 2π, basis)
fig = Figure();
showphasetight(phase .* conf2.mask, fig)
fig

# For this aberrated phase the PSF is larger

p2 = psf(conf2.ap, phase)
showarray(p2, :grays)

# More details are visible in the logarithmic scale

showarray(PhaseRetrieval.logrescale(p2))

# The `SimConfig`type is callable and, if applied to an array of proper dimensions,
# generates a psf

p2 = conf2(phase)
showarray(p2, :grays)

# ## High NA
camHR = CameraChip(;
    pixelsize=0.05um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8
)
lensHR = ImagingLens(15mm, 25mm)
imsHR = ImagingSensor(; lens=lensHR, cam=camHR)
PhaseRetrieval.numericalaperture(imsHR)

# Now we can save all these definitions in a simulation config [`SimConfig`](@ref). We also specify the wavelength here:

confHR = SimConfig("High NA", imsHR, 633nm)

# Scalar PSF
pscalar = psf(confHR.ap)
heatmap(rotr90(pscalar[503:523, 631:651]); axis=(aspect=DataAspect(),))

# Compare it with a PSF with linear polarisation in x direction.
# It should be narrow in y.
polmod = PhaseRetrieval.get_polarization_magnitudes(confHR);

pvector = PhaseRetrieval.incoherent_psf(confHR, [p for p in polmod[[:exx, :eyx, :ezx]]])
heatmap(rotr90(pvector[503:523, 631:651]); axis=(aspect=DataAspect(),))

# While we see the elliptical structure, it's not narrower in y, as in Mansuripur's example.
