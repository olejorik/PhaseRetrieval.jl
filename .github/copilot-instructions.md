# PhaseRetrieval.jl — Copilot Context

## What this package is

`PhaseRetrieval` implements the PSF/OTF forward model and phase retrieval algorithms (Gerchberg-Saxton, Alternating Projections, DRAP) for scalar and vectorial optical systems. It is the algorithmic core dependency of `Feedback14AMI`.

## Key types

| Type | Purpose |
|---|---|
| `ImagingLens` | Lens model: focal length, NA, wavelength |
| `CameraChip` | Detector: pixel size, bit depth, ROI |
| `ImagingSensor` | Combined lens + camera system |
| `SimConfig` | Forward model simulation configuration |
| `FixedExposure` | Exposure model for image simulation |
| `PRproblem` / `PDPRproblem` | Phase retrieval problem (single / phase-diversity) |
| `PRproblemSat` | Problem with saturation model |
| `GS` / `GSparam` | Gerchberg-Saxton algorithm + parameters |
| `AP` / `APparam` | Alternating Projections + parameters |
| `DRAP` / `DRAPparam` | Douglas-Rachford AP variant |
| `SHdiversity` | Shack-Hartmann diversity phase |

## Key functions

- `psf` / `subpsf` — compute PSF from pupil function (scalar or vectorial)
- `solve` — run phase retrieval (dispatches on problem + algorithm type)
- `removepiston` / `removetiptilt` — remove low-order aberrations from result
- `errhist` / `disthist` — convergence history accessors
- `gaussian_apodization` — smooth pupil edge apodization

## Vectorial PSF

`VectorialPSF.jl` implements the Richards-Wolf integral for high-NA systems with polarization handling. Used in `Feedback14AMI` via `PaxmanEstimateVectorial*`.

## Relationships

- Depends on: `PhaseBases`, `PhaseUtils`, `AlternatingProjections`, `SampledDomains`, `FFTW`
- Used by: `Feedback14AMI`

## Source layout

```
src/
    PhaseRetrieval.jl   ← module entry, exports
    hardware.jl         ← ImagingLens, CameraChip, ImagingSensor
    types.jl            ← SimConfig, exposure models
    PSF.jl              ← psf, scalar forward model
    VectorialPSF.jl     ← Richards-Wolf vectorial PSF
    PRproblem.jl        ← problem types and solve dispatch
    methods.jl          ← GS, AP, DRAP algorithms
    forwardmodel.jl     ← image simulation utilities
ShackHartmann/          ← SH sensor submodule
```
