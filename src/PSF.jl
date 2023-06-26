# Different methods of generating PSF

using FFTW

field(amplitude, phase) = amplitude .* exp.( 1im * phase)

"""
psf(amplitude, phase) -> psfimage
psf(pupilfield) -> psfimage


Calculate psf for given amplitude and phase.

Examples
====
```jldoctest

ap,_ = PhaseRetrieval.aperture(-1:.2:1,-1:.2:1,.8);

psfimage = psf(ap)

psfimage[1]

# output

17.720305330423212

```
"""
psf(field::AbstractArray) = abs.(toimageplane(field)) .^2

psf(amplitude, phase) = psf(field(amplitude, phase))

subpsf(field::AbstractArray,Q::Integer) = psf(subdivide_sum(field,Q))
subpsf(amplitude, phase,Q::Integer) = psf(subdivide_sum(field(amplitude, phase),Q))


toimageplane(field) = fftshift(fft(ifftshift(field)))



