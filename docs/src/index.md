```@meta
CurrentModule = PhaseRetrieval
DocTestSetup = quote
    using PhaseRetrieval
end
```


# PhaseRetrieval.jl Documentation

```@contents
```

## About the package
This package is devoted to the forward and inverse problems of Phase Retrieval (PR).

## Basic usage
### Forward model 
Let's set up a simulation environment matching the following hardware set up: a beam with a footprint of 10 mm diameter is focused with a lens of 300 mm focal length and the PSF is registered with UI-1540 camera, which we combine in one structure called `ImagingSensor`(@ref PhaseRetrieval.ImagingSensor)

```@example psf
using PhaseRetrieval
lens = PhaseRetrieval.ImagingLens(300mm, 25mm)
cam = PhaseRetrieval.CameraChip(pixelsize = 5.2um, imagesize = (1280, 1024), bitdepth = 8, channelbitdepth = 8)
ims = PhaseRetrieval.ImagingSensor(lens = lens, cam = cam)
```

Now we can save all these defintitions in a simulatin configg
```@example psf
conf1 = SimConfig("full_aperture", ims, 633nm)
```

This creates aperture array of correct dimensions which is suitable for generation of a PSF using Fourier methods. If we check near the central pixel, we'll see that for this configuration the PSF is almost one pxels wide
```@example psf
p = psf(conf1.ap)
using CairoMakie # hide
CairoMakie.activate!(type = "svg") # hide
heatmap(rotr90(p[503:523,631:651]), axis = (aspect = DataAspect(), ))
```
Indeed, the Airy pattern should be about 9 microns wide
```@example
print("Airy size is 1.22λ/NA =",  1.22*632nm *300mm/25mm /um, " μm")
```

We might thus want to consider a smaller numerical aperture:

```@example psf
lens2= PhaseRetrieval.diaphragm(lens, 10mm)
ims2 = PhaseRetrieval.ImagingSensor(lens = lens2, cam = cam)
conf2 = SimConfig("10mm aperture", ims2, 633nm)
p2 = psf(conf2.ap)
heatmap(rotr90(p2[503:523,631:651]), axis = (aspect = DataAspect(), ))
```

## Types and Functions

```@autodocs
Modules = [PhaseRetrieval]
```

## Index
```@index
```
