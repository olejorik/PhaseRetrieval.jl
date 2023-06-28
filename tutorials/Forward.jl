# ```@meta
# CurrentModule = PhaseRetrieval
# DocTestSetup = quote
#     using PhaseRetrieval
# end
# ```

# ### Forward model 

# For the PR problem, forward model is simulation of a realistic readout of a camera of a PSF (or an extended object) under some predefined conditions.

# Let's set up a simulation environment matching the following hardware set up: a beam with a footprint of 1 inch (25 mm) diameter is focused with a lens of 300 mm focal length and the PSF is registered with [UI-1240 camera](https://en.ids-imaging.com/store/products/cameras/ui-1240le.html).
# For both lens and camera, we can use structures with self-explanatory names
# [`ImagingLens`](@ref) and [`CameraChip`](@ref), which we combine in one structure called [`ImagingSensor`](@ref):

using PhaseRetrieval
lens = PhaseRetrieval.ImagingLens(300mm, 25mm)
cam = PhaseRetrieval.CameraChip(;
    pixelsize=5.3um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8
)
ims = PhaseRetrieval.ImagingSensor(; lens=lens, cam=cam)

# Now we can save all these definitions in a simulation config [`SimConfig`](@ref). We also specify the wavelength here:

conf1 = SimConfig("full_aperture", ims, 633nm)

# This creates aperture array of correct dimensions which is suitable for generation of a PSF using Fourier methods. If we check near the central pixel, we'll see that for this configuration the PSF is almost one pixel wide

p = psf(conf1.ap)
using CairoMakie # hide
CairoMakie.activate!(; type="png") # hide
heatmap(rotr90(p[503:523, 631:651]); axis=(aspect=DataAspect(),))

# Indeed, the Airy pattern should be about 9 microns wide

print("Airy size is 1.22λ/NA = ", 1.22 * 632nm * conf1.f / conf1.d / um, " μm")

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

ims = ImagingSensor(; cam=cam = camerasdict["UI1240"], lens=lensesdict["F300A25"])

# #### `SimConfig`
# `SimConfig` contains necessary information for simulations.

fieldnames(typeof(conf2))

# For instance, it contains the aperture mask.

using PhasePlots
showarray(conf2.ap)

# The dimensions of the mask correspond to the dimensions of the sampled image plane, but the overall size corresponds to the inverse of the pixel size. This information is contained in `dualroi` field and can be used to construct the Zernike basis.

conf2.dualroi

using PhaseBases
basis = ZernikeBW(conf2.dualroi, conf2.d, 10);
showphase(basis.elements[15] .* conf2.mask)
current_figure()

# Or the same picture without unnecssary information (by default all phases will be shown scaled to (-π. π])

showphasetight(basis.elements[15] .* conf2.mask)
current_figure()

# This is a combination of some low-order Zernike polynomials

phase = compose(basis, [4, 6, 15, 16], [2, 1, 0.4, 0.3] * 2π)
fig = Figure();
showphasetight(phase .* conf2.mask, fig)
fig

# For this aberrated phase the PSF is larger

p2 = psf(conf2.ap, phase)
showarray(p2, :grays)

# More details are visible in the logarithmic scale

showarray(PhaseRetrieval.logrescale(p))

# The `SimConfig`type is callable and, if applied to an array of a proper dimensions,
# generates a psf

p2 = conf2(phase)
showarray(p2, :grays)

# # Inverse problem
# The goal of the inverse problem is from given PSF and `SimConfig`  
# to restore the unknow phase.

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
# ![Output after 1500 iterations](assets/PR_DRAP1500.png)

#  The phase is perfectly reconstructed now, and 
#  so here is the main problem of the AP-based PR algorthms --- they require quite
#  a long time to converge, even for the noiseless data.

# Let's try to have a smaller crop
center = [size(p2)...].÷2
crophw=64
pcrop = p2[CartesianIndex(center...) - CartesianIndex(crophw, crophw):CartesianIndex(center...) + CartesianIndex(crophw-1, crophw-1)]
showarray(pcrop)

# Construct corresponding sim config and see how it works
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

# You can try to change slightly the values of `β` above and see that algortihm
# might converge to another solution. This is another problem of AP-based algorithms.