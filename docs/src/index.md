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
Let's set up a simulation environment matching the following hardware set up: a beam with a footprint of 10 mm diameter is focused with a lens of 300 mm focal length and the PSF is registered with UI-1540 camera.

```@example
using PhaseRetrieval
lens = PhaseRetrieval.ImagingLens(300mm, 25mm)
```

## Types and Functions

```@autodocs
Modules = [PhaseRetrieval]
```

## Index
```@index
```
