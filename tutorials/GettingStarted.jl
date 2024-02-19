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
cam = CameraChip(; pixelsize=5.3um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8)
ims = ImagingSensor(; lens=lens, cam=cam)

# Now we can save all these definitions in a simulation config [`SimConfig`](@ref).
# We also specify the wavelength here:

conf1 = SimConfig("full_aperture", ims, 633nm)

# This creates a simulation configuration for the generation of a PSF using [`Fourier`](@ref) methods.
# If we check near the central pixel, we'll see that for this configuration the PSF is
# almost one pixel wide.
# To crop central part of the image, we'll use `crop` function from `PhaseUtils` package.

p = psf(conf1)
using CairoMakie # hide
CairoMakie.activate!(; type="png") # hide
using PhaseUtils
heatmap(rotr90(crop(p, 9)); axis=(aspect=DataAspect(),)) # show central 9 pixels

# Indeed, the Airy pattern should be about 9 microns wide

print("Airy size is 1.22λ/NA = $(airysize(conf1)/ um) μm")

# We might thus want to consider a smaller numerical aperture:

lens2 = PhaseRetrieval.diaphragm(lens, 10mm)
ims2 = PhaseRetrieval.ImagingSensor(; lens=lens2, cam=cam)
conf2 = SimConfig("10mm aperture", ims2, 633nm)
p2 = psf(conf2)
heatmap(rotr90(crop(p2, 9)); axis=(aspect=DataAspect(),))

# ### Faster creation of an `ImagingSensor`
# Some often used cameras are saved in [`camerasdict`](@ref) and [`lensesdict`](@ref) dictionaries

keys(camerasdict)

# So the imaging sensor can be created as

ims = ImagingSensor(; cam=camerasdict["UI1240"], lens=lensesdict["F300A25"])

# ### Quantisation and exposure level
# By default, the returned PSF approximates the output of a camera with
# finite bit resolution (8 bits in our case):

eltype(p2)

# One can use `Images.jl` package to obtain images as saved by the camera:
# The whole frame:
using ImageCore
save("psf2.png", Gray.(p2))

# ![](psf2.png)

# Crop of the central part
Gray.(crop(p2, 21))

# By default, the returned PSF is called between 0 and 1 ([`AutoExposure`](@ref) feature).
# This can be changed by passing additional parameters to [`psf`](@ref) function.
# Here is an example of a PSF with 4 times longer exposure:
psf2_sat4 = psf(conf2; exposure=AutoExposure(4))
Gray.(crop(psf2_sat4, 21))

# The same psf in logarithmic scale:
Gray.(logrescale(float.(crop(psf2_sat4, 21))))

# And without quantisation:
psf2_sat4_float = psf(conf2; exposure=AutoExposure(4), quantize=false)
Gray.(logrescale(crop(psf2_sat4_float, 21)))

# ### Adding the phase aberration
# Now we can add some phase to our configuration.
# To add a modal phase represented by Zernike polynomials, we need to
# create the basis first. This creates it for the first 10 radial orders:

using PhaseBases
z10 = ZernikeBW(conf2, 10);

# The basis now contains 60 Zernike polynomials (in Born and Wolf form,
# normed by *rms* value),
# numbered in OSA/ANSI indexes. Elements of basis can be accessed as follows:
# - by double indexing:

using PhasePlots
showphasetight(z10(; m=10, n=10) .* conf2.mask)[1]

# - by single index:
showphasetight(z10(12) .* conf2.mask)[1]

# This is a combination of some low-order Zernike polynomials
phase = ModalPhase([4, 6, 15, 16], [2, 1, 0.4, 0.3] * 2π, z10)
showphasetight(phase .* conf2.mask)[1]

# To apply this phase to the simulation configuration, apply it literally:
phase(conf2)

# and check the PSF
p2 = psf(conf2)
save("psf2_ab.png", Gray.(p2));
nothing #hide
# ![](psf2_ab.png)

# For this aberrated phase the PSF is larger and
# more details are visible in the logarithmic scale

showarray(logrescale(float(p2)))

# The `SimConfig`type is callable and, if applied to an array of proper dimensions,
# generates a psf

p2 = conf2(phase)
showarray(p2, :grays)

# ### Adding some diversity
# This creates a defocus of 1λ/4 amplitude:
defocus = 2π / 4 * z10(; n=2, m=0)
showphase(collect(defocus .* conf2.mask))[1]

#  Or we could better create diversities corressponding in the depth-of-focus lentgths,
# which is approximately the same
defocus = ZonalPhase(throughfocus(conf2, doflength(conf2)))
showphase(collect(-(defocus .+ π / 2) .* conf2.mask))[1]

# And now we add several diversities to our `conf2`
for k in [-2, -1, 1, 2]
    conf2.diversity[string(k)] = P = k * defocus
end;

div_psf = diversed_psfs(conf2)
for (p, d) in zip(div_psf, vcat(["0"], collect(keys(conf2.diversity))))
    fig, ax, hm = showarray(logrescale(float(crop(p, 256))))
    ax.title = "PSF #" * d
    save("psf2_div_fig" * d * ".png", fig)
end
# ![](psf2_div_fig-2.png)
# ![](psf2_div_fig-1.png)
# ![](psf2_div_fig0.png)
# ![](psf2_div_fig1.png)
# ![](psf2_div_fig2.png)

# ## [Inverse problem](@id tutorial-inverse)
# The goal of the inverse problem is from the given PSF and `SimConfig`
# to restore the unknown phase.

# In mathematics, the (2D) phase retrieval (PR) problem can be described as
# ```math
# \text{find }x \mathbb{C}^{M \times N} \text{s.t.}\\
# \abs{x} = a, \abs{\F x} = A
# ```
# for some real arrays $a$ and $A.$
#
# As for the PSF $p$ one has $p = \abs{\F (a e^{i \phi})}$, we can from the PRproblem by passing to `PRproblem` object two arrays corresponding to the square root of the intensities in the pupil and focal planes.
# This PR problem we try to solve using successive applications of the DRAP algorithm with β = 0.9 and Alternating projections algorithm with the default number of iterations.

using FFTW
a = sqrt.(Float64.(conf2.ap))
A = sqrt.(collect(Float64, p2))
pr = PRproblem(a, A)
sol = solve(pr, (DRAPparam(; β=0.9, keephistory=true), APparam(; maxϵ=0.001)))
showphasetight(fftshift(angle.(sol[1])) .* conf2.mask)[1]

# This doesn't work at the moment.
#  It might work if we let it iterate further
# ```julia
# sol = solve(pr, (DRAPparam(β = 0.9,keephistory = true, maxit =1500), APparam(maxϵ = 0.001)))
# showphasetight(fftshift(angle.(sol[1])) .* conf2.mask)[1]
# ```
# ![Output after 1500 iterations](../assets/PR_DRAP1500.png)

#  The phase is perfectly reconstructed now, and so here is the main problem of the AP-based PR algorithms --- they require quite a long time to converge, even for the noiseless data.

# Here is an attempt with the initial guess given by the simple subset method
th = 0.6 * maximum(A)
xth = copy(A)
xth[A .> th] .= 1
xth[A .<= th] .= 0
showarray(xth)

alg = DRAPparam(; x⁰=ifft(xth), β=0.95, keephistory=true, maxit=200)
sol = solve(pr, (alg, APparam(; maxit=50)))
showphasetight(fftshift(angle.(sol[1])) .* conf2.mask)[1]

# Let's try to have a smaller crop
cropw = 128
pcrop = crop(p2, cropw)
showarray(pcrop)

# Construct the corresponding sim config and see how it works
# TODO wrap all this in functions
ims2crop = PhaseRetrieval.ImagingSensor(; lens=lens2, cam=PhaseRetrieval.roi(cam, cropw))
conf2crop = SimConfig("10mm aperture", ims2crop, 633nm)
a = sqrt.(Float64.(conf2crop.ap))
A = sqrt.(collect(Float64, pcrop))
pr = PRproblem(a, A)
sol = solve(pr, (DRAPparam(; β=0.9, keephistory=true, maxit=550), APparam(; maxit=10)))
# sol = solve(pr, (DRAPparam(β = 0.9,keephistory = true), APparam(maxϵ = 0.001)))
showphasetight(fftshift(angle.(sol[1])) .* conf2crop.mask)[1]

# You can try to change slightly the values of `β` above and see that algorithm
# might converge to another solution. This is another problem of AP-based algorithms.
sol = solve(pr, (DRAPparam(; β=0.91, keephistory=true, maxit=500), APparam(; maxit=10)))
showphasetight(fftshift(angle.(sol[1])) .* conf2crop.mask)[1]

# Douglas-Rachford is known to eventually find the solution if you run it long enough:
# For instance, starting with `b`=0.91 would require about 20000 iterations to converge:
sol = solve(pr, (DRAPparam(; β=0.91, keephistory=true, maxit=20000), APparam(; maxit=100)))
showphasetight(fftshift(angle.(sol[1])) .* conf2crop.mask)[1]

# Fortunately, julia is fast, so the calculations of 20K iterations take less than a minute.

##
# ## Using phase diversity
#  To use phase-diverse Phase Retrieval, we need to construct the phase diversities for the cropped configuration
z10 = ZernikeBW(conf2crop, 10);
# Then we know, that cropped PSFs correspond to different values of the same defocus
## defocus = 2π / 4 * z10(; n=2, m=0)
defocus = ZonalPhase(throughfocus(conf2crop, doflength(conf2crop)))
phases = [collect(k * defocus) for k in [0, 1, 2, -2, -1]] #TODO this should be automatized
# Crop the corresponding PSFs
div_psf_crop = crop.(div_psf, cropw)
save("psf_crop.png", Gray.(mosaicview(div_psf_crop; nrow=1, npad=5, fillvalue=1)))

# ![](psf_crop.png)

# Now the set corresponding to the pupil phase will be phase-diversed set

a = sqrt.(Float64.(conf2crop.ap))
A = [sqrt.(collect(Float64, p)) for p in div_psf_crop]
pr = PDPRproblem(a, A, phases)


# As expected, using phase diversities accelerates the phase retrieval a lot
sol = solve(pr, (DRAPparam(; β=0.9, keephistory=true, maxit=50), APparam(; maxit=50)))
showphasetight(fftshift(angle.(sol[1][:, :, 1])) .* conf2crop.mask)[1]

# And now it should be also less sensitive to the parameter choice.
# But they do affect the speed of the convergence
sol = solve(pr, (DRAPparam(; β=0.5, keephistory=true, maxit=50), APparam(; maxit=50)))
showphasetight(fftshift(angle.(sol[1][:, :, 1])) .* conf2crop.mask)[1]

# Try to unwrap the phase
ph = fftshift(angle.(sol[1][:, :, 1]))
ph_u = unwrap_LS(ph, conf2crop.ap)
fig, ax, hm = showarray(bboxview(ph_u .* conf2crop.mask))
Colorbar(fig[1, 2], hm)
fig

# Now the phase is unwrapped and can be decomposed by the basis functions
restored_coef = PhaseBases.decompose(ph_u, z10)
#  and compare the restored coefficients with the original ones
fig = Figure();
ax = Axis(fig[1, 1]);
scatter!(phase.coef; label="original")
scatter!(restored_coef; label="restored")
ax.title = "Original vs restored Zernike coefficients"
axislegend()

err_coef = restored_coef .- phase.coef
err_coef[1] = 0
ax, _ = scatter(fig[2, 1], err_coef)
ax.title = "Coefficient restoration error"
fig

# The error is quite small, although we see it increasing for high-order aberrations.
#  We can compare the phases themselves and see the error appearing on the edges.
ph_u .-= restored_coef[1]
phase_crop = ModalPhase(phase.coef, z10)
err = (ph_u .- phase_crop)
fig, ax, hm = showarray(bboxview(err .* conf2crop.mask))
Colorbar(fig[1, 2], hm)
ax.title = "Phase restoration error, rms = $(round(maskedrmse(err, z10.ap), digits = 3))"
fig

# Try with high-res phase
z10 = ZernikeBW(conf2, 10);
defocus = 2π / 4 * z10(; n=2, m=0)
defocus = ZonalPhase(throughfocus(conf2, doflength(conf2)))
phases = [collect(k * defocus) for k in [0, 1, 2, -2, -1]] #TODO this should be automatized

a = sqrt.(Float64.(conf2.ap))
A = [sqrt.(collect(Float64, p)) for p in div_psf]

pr = PDPRproblem(a, A, phases)
sol = solve(pr, (DRAPparam(; β=0.5, keephistory=true, maxit=50), APparam(; maxit=50)))
showphasetight(fftshift(angle.(sol[1][:, :, 1])) .* conf2.mask)[1]

# We again unwrap the phase
ph = fftshift(angle.(sol[1][:, :, 1]))
ph_u = unwrap_LS(ph, conf2.ap)
fig, ax, hm = showarray(bboxview(ph_u .* conf2.mask))
ax.aspect = 1
Colorbar(fig[1, 2], hm)
fig

# Now the phase is unwrapped and can be decomposed by the basis functions
restored_coef = PhaseBases.decompose(ph_u, z10)
#  and compare the restored coefficients with the original ones
fig = Figure();
ax = Axis(fig[1, 1]);
scatter!(phase.coef; label="original")
scatter!(restored_coef; label="restored")
ax.title = "Original vs restored Zernike coefficients"
axislegend()

err_coef = restored_coef .- phase.coef
err_coef[1] = 0
ax, _ = scatter(fig[2, 1], err_coef)
ax.title = "Coefficient restoration error"
fig
# The error is error is about the same.

#  We can compare the phases themselves
ph_u .-= restored_coef[1]
err = (ph_u .- phase)
fig, ax, hm = showarray(bboxview(err .* conf2.mask));
ax.aspect = 1
Colorbar(fig[1, 2], hm)
ax.title = "Phase rmse = $(PhaseUtils.maskedphasermse(ph_u, phase, conf2.ap)/(2π)) λ"
fig

# The quite big error in the coefficient restoration compared with the error magnitude in the phases _is not_ the numerical error of decomposition by Zernikes.
# Here is, for instance, the result of the decomposition of the input phase:
rest = PhaseBases.decompose(collect(phase), z10)
scatter(phase.coef)
scatter!(rest)
current_figure()

# But the error grows with the index of the Zernike:
scatter(rest .- phase.coef; axis=(title="RMS error = $( norm(rest .- phase.coef))",))
