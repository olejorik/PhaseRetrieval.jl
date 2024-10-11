# Different methods used for PSF simulations
# Concrete methods are intentionally not exported, so we can discriminate for `Fourier` method for PSF simulation and for other types of problems.
# Abstract method is exported for the interface definitions.

module PSFMethods

export PSFMethod

"""
    PSFMethod

Abstract method for simulating PSF.
"""
abstract type PSFMethod end

"""
    PSFMethods.Fourier()

FFT-based method of simulation of a PSF.

See also `psf`, `SimConfig`.
"""
struct Fourier <: PSFMethod end
struct ChirpZ <: PSFMethod end
struct MVM <: PSFMethod end

end


using .PSFMethods

#  defaults and interfaces
const default_forward_method = PSFMethods.Fourier()


"""
    toimageplane(field [, method = PSFMethods.Fourier()])

Calculate the field in the focal plane from the field in the pupil plane using `method`.

Available methods:
 - `PSFMethods.Fourier()`
 - `ChirpZ`  -- not implemented
 - `MVM` -- not implemented
"""
toimageplane(field) = toimageplane(field, default_forward_method)

(m::PSFMethod)(f::AbstractField) = m(f, polarisation(f))
(m::PSFMethod)(f::AbstractField, p::AbstractPolarization) = error(
    "$(typeof(m)) method for field $(typeof(f)) and polarisation $(typeof(p)) is not implemented",
)

# Scalar polarization -- just calculate the field in the focal plane
(m::PSFMethod)(f::AbstractField, ::ScalarPolarization) =
    toimageplane(collect(complex_amplitude(f)), m)

# Linear polarization -- calculate three field components for x- and y-polarization components and sum them coherently. We assume that the poalrization modulations are given as named tuple
function (m::PSFMethod)(f::AbstractField, p::LinearPolarization, pm)
    s, c = sincos(p.orientation_angle)
    xcomp = toimageplane(collect(complex_amplitude(f) .* (c * pm.exx + s * pm.exy)), m)
    ycomp = toimageplane(collect(complex_amplitude(f) .* (c * pm.eyx + s * pm.eyy)), m)
    zcomp = toimageplane(collect(complex_amplitude(f) .* (c * pm.ezx + s * pm.ezy)), m)

    return (xcomp, ycomp, zcomp)

end

# Fourier defined for an array

toimageplane(field::AbstractArray, ::PSFMethods.Fourier) = fftshift(fft(ifftshift(field)))
