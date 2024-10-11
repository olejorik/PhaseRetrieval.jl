# Collection of low-level computational algorithms specific for PSF simulation and phase retrieval
#


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
psf(field::AbstractArray) = myabs2.(toimageplane(field))
# psf(field::AbstractArray) = abs.(toimageplane(field)) .^2

psf(amplitude, phase) = psf(field(amplitude, phase))

subpsf(field::AbstractArray, Q::Integer) = psf(subdivide_sum(field, Q))
subpsf(amplitude, phase, Q::Integer) = subpsf(field(amplitude, phase), Q)

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
intensity(field) = mappedarray(myabs2, field)
