# ```@meta
# CurrentModule = PhaseRetrieval
# DocTestSetup = quote
#     using PhaseRetrieval
# end
# ```

# # Getting started
# The goal of this package is to provide tools for the forward and inverse problems of the 
# Phase Retrieval (PR).
# The examples below show how to set up a simulation environment and how to form a PR problem from 
# given measured PSF and an aperture.



# ```@contents
# Pages = [
#     "GettingStarted.md",
# ]
# Depth = 2
# ```

# ## [Forward model](@id tutorial-forward) 

# For the PR problem, the forward model is a simulation of a 
# realistic readout of a camera of a PSF (or an extended object) under some predefined conditions.

# ### Initialising the simulation
# Let's set up a simulation environment matching the following hardware setup: 
# a beam with a footprint of 1 inch (25 mm) diameter is focused with a lens of 300 mm 
# focal length, and the PSF is registered with [UI-1240 camera](https://en.ids-imaging.com/store/products/cameras/ui-1240le.html).
# For both lens and camera, we can use structures with self-explanatory names
# [`ImagingLens`](@ref) and [`CameraChip`](@ref), which we combine in one structure 
# called [`ImagingSensor`](@ref):

using PhaseRetrieval
lens = ImagingLens(300mm, 25mm)
cam = CameraChip(;
    pixelsize=5.3um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8
)
ims = ImagingSensor(lens=lens, cam=cam)

# Now we can save all these definitions in a simulation config [`SimConfig`](@ref). 
# We also specify the wavelength here:

conf1 = SimConfig("full_aperture", ims, 633nm)

# This creates a simulation configuration for the generation of a PSF using [`Fourier`](@ref) methods. 
# If we check near the central pixel, we'll see that for this configuration the PSF is 
# almost one pixel wide

p = psf(conf1)
using CairoMakie # hide
CairoMakie.activate!(; type="png") # hide
heatmap(rotr90(p[503:523, 631:651]); axis=(aspect=DataAspect(),))

# Indeed, the Airy pattern should be about 9 microns wide

print("Airy size is 1.22λ/NA = $(airysize(conf1)/ um) μm")

# We might thus want to consider a smaller numerical aperture:

lens2 = PhaseRetrieval.diaphragm(lens, 10mm)
ims2 = PhaseRetrieval.ImagingSensor(; lens=lens2, cam=cam)
conf2 = SimConfig("10mm aperture", ims2, 633nm)
p2 = psf(conf2)
heatmap(rotr90(p2[503:523, 631:651]); axis=(aspect=DataAspect(),))

# ### Faster creation of an `ImagingSensor`
# Some often used cameras are saved in [`camerasdict`](@ref) and [`lensesdict`](@ref) dictionaries

keys(camerasdict)

# So the imaging sensor can be created as 

ims = ImagingSensor(; cam = camerasdict["UI1240"], lens=lensesdict["F300A25"])

# ### Quantisation and exposure level
# By default, the returned PSF approximates the output of a camera with 
# finite bit resolution (8 bits in our case):

eltype(p2)

# One can use `Images.jl` package to obtain images as saved by the camera:
# The whole frame:
using ImageCore
p2 .|> Gray

# Crop of the central part
p2[503:523, 631:651] .|> Gray

# By default, the returned PSF is called between 0 and 1 ([`AutoExposure`](@ref) feature).
# This can be changed by passing additional parameters to [`psf`](@ref) function.
# Here is an example of a psf with 4 times longer exposure:
psf2_sat4 = psf(conf2, exposure =  AutoExposure(4))
psf2_sat4[503:523, 631:651]   .|> Gray

# The same psf in logarithmic scale:
psf2_sat4[503:523, 631:651] .|> float |> logrescale  .|> Gray

# And ithout quantisation:
psf2_sat4_float = psf(conf2, exposure =  AutoExposure(4), quantize = false)
psf2_sat4_float[503:523, 631:651] |> logrescale  .|> Gray


# ### Adding the phase aberration
# Now we can add some phase to out configuration.
# To add a modal phase represented by Zernike polynomials, we nee to
# create the basis first. This creates it for the first 10 radial orders:

using PhaseBases
z10 = ZernikeBW(conf2, 10);


# The basis now contains 60 Zernike polynomials (in Born and Wolf form,
# normed by *rms* value),
# numbered in OSA/ANSI indexes. Elements of basis can be accessed as follows:
# - by double indexing:

using PhasePlots
showphasetight(z10(m=10,n=10) .* conf2.mask)[1]

# - by single index:
showphasetight(z10(12) .* conf2.mask)[1]


# This is a combination of some low-order Zernike polynomials
phase = ModalPhase([4, 6, 15, 16], [2, 1, 0.4, 0.3] * 2π, z10)
showphasetight(phase .* conf2.mask)[1]

# To apply this phase to the simulation configuration, apply it literally:
phase(conf2)

# and chenck the PSF
p2 = psf(conf2)
p2 .|> Gray

# For this aberrated phase the PSF is larger and 
# more details are visible in the logarithmic scale

showarray(logrescale(float(p2)))

# The `SimConfig`type is callable and, if applied to an array of proper dimensions,
# generates a psf

p2 = conf2(phase)
showarray(p2, :grays)

# ## [Inverse problem](@id tutorial-inverse)
# The goal of the inverse problem is from the given PSF and `SimConfig`  
# to restore the unknown phase.

# As at this stage the problem is already reduced to its numerical equivalent 
# ```math
# ```

using FFTW
using AlternatingProjections
# pr = PRproblem(conf2.ap, p)
a = fftshift(sqrt.(Float64.(conf2.ap)))
A = fftshift(sqrt.(collect(Float64, p2))) 
N = sqrt(sum(abs2, A))
n = sqrt(sum(abs2, a))
A = A ./ N .* n
pr = TwoSetsFP(
                ConstrainedByAmplitude(a), FourierTransformedSet(ConstrainedByShape(A))
            )
sol = solve(pr, (DRAPparam(β = 0.9,keephistory = true, ), APparam(maxϵ = 0.001)))
showphasetight(fftshift(angle.(sol[1])) .* conf2.mask); current_figure()

# This doesn't work at the moment. 
#  It might work if we let it iterate further
# ```julia
# sol = solve(pr, (DRAPparam(β = 0.9,keephistory = true, maxit =1500), APparam(maxϵ = 0.001)))
# showphasetight(fftshift(angle.(sol[1])) .* conf2.mask); current_figure()
# ```
# ![Output after 1500 iterations](../assets/PR_DRAP1500.png)

#  The phase is perfectly reconstructed now, and 
#  so here is the main problem of the AP-based PR algorithms --- they require quite
#  a long time to converge, even for the noiseless data.

# Let's try to have a smaller crop
center = [size(p2)...].÷2
crophw=64
pcrop = p2[CartesianIndex(center...) - CartesianIndex(crophw, crophw):CartesianIndex(center...) + CartesianIndex(crophw-1, crophw-1)]
showarray(pcrop)

# Construct the corresponding sim config and see how it works
# TODO wrap all this in functions
ims2crop = PhaseRetrieval.ImagingSensor(; lens=lens2, cam=PhaseRetrieval.roi(cam, 2crophw))
conf2crop = SimConfig("10mm aperture", ims2crop, 633nm)
a = fftshift(sqrt.(Float64.(conf2crop.ap)))
A = fftshift(sqrt.(collect(Float64, pcrop ))) 
N = sqrt(sum(abs2, A))
n = sqrt(sum(abs2, a))
A = A ./ N .* n
pr = TwoSetsFP(
                ConstrainedByAmplitude(a), FourierTransformedSet(ConstrainedByShape(A))
            )
sol = solve(pr, (DRAPparam(β = 0.9,keephistory = true, maxit=450), APparam(maxit = 10)))
# sol = solve(pr, (DRAPparam(β = 0.9,keephistory = true), APparam(maxϵ = 0.001)))
showphasetight(fftshift(angle.(sol[1])) .* conf2crop.mask); current_figure()

# You can try to change slightly the values of `β` above and see that algorithm
# might converge to another solution. This is another problem of AP-based algorithms.