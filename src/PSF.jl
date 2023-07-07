# Different methods of generating PSF

using FFTW
using MappedArrays

abstract type PSFmethod end

"""
    Fourier()

FFT-based method of simulation of a PSF.
"""
struct Fourier <: PSFmethod end
struct ChirpZ <: PSFmethod end
struct MVM <: PSFmethod end

abstract type PSFExposure end
struct AutoExposure <:PSFExposure
    scale::Float64
end

AutoExposure() = AutoExposure(1)

export AutoExposure

field(amplitude, phase) = amplitude .* exp.( 1im * phase)

"""
    psf(amplitude, phase) -> psfimage
    psf(pupilfield) -> psfimage


Calculate psf for given amplitude and phase. Both PSF and pupil arrays are assumed
    to have the origin in their center pixel (`(N+1)÷2` for dimension size `N`).

Examples
====
```jldoctest
julia> ap,_ = PhaseRetrieval.aperture(-1:.2:1,-1:.2:1,.8); ap
11×11 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> Int.(round.(psf(ap)))
11×11 Matrix{Int64}:
 18  3   1   2    0    1    0   2   1  3  18
  3  0   6   5    0    1    0   5   6  0   3
  1  6  11   4    1    5    1   4  11  6   1
  2  5   4   1   20   38   20   1   4  5   2
  0  0   1  20   78  118   78  20   1  0   0
  1  1   5  38  118  169  118  38   5  1   1
  0  0   1  20   78  118   78  20   1  0   0
  2  5   4   1   20   38   20   1   4  5   2
  1  6  11   4    1    5    1   4  11  6   1
  3  0   6   5    0    1    0   5   6  0   3
 18  3   1   2    0    1    0   2   1  3  18

```

```jldoctest
julia> psf(ones(3,3),zeros(3,3))
3×3 Matrix{Float64}:
 0.0   0.0  0.0
 0.0  81.0  0.0
 0.0   0.0  0.0
 
julia> psf(ones(Complex{Float64},4,4))
4×4 Matrix{Float64}:
 0.0  0.0    0.0  0.0
 0.0  0.0    0.0  0.0
 0.0  0.0  256.0  0.0
 0.0  0.0    0.0  0.0

```
"""
# psf(field::AbstractArray) = abs2.(toimageplane(field)) 
psf(field::AbstractArray) = abs.(toimageplane(field)) .^2

psf(amplitude, phase) = psf(field(amplitude, phase))

subpsf(field::AbstractArray,Q::Integer) = psf(subdivide_sum(field,Q))
subpsf(amplitude, phase,Q::Integer) = psf(subdivide_sum(field(amplitude, phase),Q))


"""
    toimageplane(field [, method = Fourier()])

Calculate the field in the focal plane from the field in the pupil plane using `method`.

Available methods:
 - `Fourier()`
 - `ChirpZ`  -- not implemented
 - `MVM` -- not implemented
"""
toimageplane(field, ::Fourier) = fftshift(fft(ifftshift(field)))
toimageplane(field) = toimageplane(field, Fourier()) 

"""
    intensity(field)

Calculate intensity of the field.

# Examples

```jldoctest
julia> ap,_ = PhaseRetrieval.aperture(-1:.2:1,-1:.2:1,.8);

julia> imf = PhaseRetrieval.toimageplane(ap);

julia> PhaseRetrieval.intensity(imf) .|> round
11×11 Matrix{Float64}:
 18.0  3.0   1.0   2.0    0.0    1.0    0.0   2.0   1.0  3.0  18.0
  3.0  0.0   6.0   5.0    0.0    1.0    0.0   5.0   6.0  0.0   3.0
  1.0  6.0  11.0   4.0    1.0    5.0    1.0   4.0  11.0  6.0   1.0
  2.0  5.0   4.0   1.0   20.0   38.0   20.0   1.0   4.0  5.0   2.0
  0.0  0.0   1.0  20.0   78.0  118.0   78.0  20.0   1.0  0.0   0.0
  1.0  1.0   5.0  38.0  118.0  169.0  118.0  38.0   5.0  1.0   1.0
  0.0  0.0   1.0  20.0   78.0  118.0   78.0  20.0   1.0  0.0   0.0
  2.0  5.0   4.0   1.0   20.0   38.0   20.0   1.0   4.0  5.0   2.0
  1.0  6.0  11.0   4.0    1.0    5.0    1.0   4.0  11.0  6.0   1.0
  3.0  0.0   6.0   5.0    0.0    1.0    0.0   5.0   6.0  0.0   3.0
 18.0  3.0   1.0   2.0    0.0    1.0    0.0   2.0   1.0  3.0  18.0


```

"""
intensity(field) = mappedarray(abs2,field)

function (cam::CameraChip)(field, exposure::PSFExposure = AutoExposure()) 
    imf = intensity(field)
    maxint = maximum(imf) / exposure.scale
    storagetype = getstoragetype(cam)
    # storagetype = N4f12 #debug -- it's three times faster comapred with the unnknown type
    return mappedarray(scaleminmax(storagetype,0,maxint), imf)
end



