# The Field structure should represent the optical field in a plane. It should contain the values of the complex amplitude, wavelength, and the polarization.
#

export ScalarComplexAmplitude,
    ScalarPolarization,
    LinearPolarization,
    RandomPolarization,
    CircularPolarization,
    EllipticalPolarization,
    SampledField,
    LaserBeam


abstract type AbstractScalarComplexAmplitude end

"""
    ScalarComplexAmplitude <: AbstractScalarComplexAmplitude

Structure containing separately in its fields amplitude and phase.

"""
struct ScalarComplexAmplitude{TA,TPH} <: AbstractScalarComplexAmplitude
    amplitude::TA
    phase::TPH
end

amplitude(ca::ScalarComplexAmplitude) = ca.amplitude
phase(ca::ScalarComplexAmplitude) = ca.phase
collect(ca::ScalarComplexAmplitude) = field(amplitude(ca), phase(ca))

struct VectorComplexAmplitude{TA,TPH} <: AbstractScalarComplexAmplitude
    amplitude_x::TA
    amplitude_y::TA
    amplitude_z::TA
    phase_x::TPH
    phase_y::TPH
    phase_z::TPH
end

amplitude(ca::VectorComplexAmplitude, ::Val{:x}) = ca.amplitude_x
amplitude(ca::VectorComplexAmplitude, ::Val{:y}) = ca.amplitude_y
amplitude(ca::VectorComplexAmplitude, ::Val{:z}) = ca.amplitude_z
amplitude(ca::VectorComplexAmplitude) = (ca.amplitude_x, ca.amplitude_y, ca.amplitude_z)

phase(ca::VectorComplexAmplitude, ::Val{:x}) = ca.phase_x
phase(ca::VectorComplexAmplitude, ::Val{:y}) = ca.phase_y
phase(ca::VectorComplexAmplitude, ::Val{:z}) = ca.phase_z
phase(ca::VectorComplexAmplitude) = (ca.phase_x, ca.phase_y, ca.phase_z)
collect(ca::VectorComplexAmplitude) = field(amplitude(ca), phase(ca))


abstract type AbstractPolarization end

"""
    LinearPolarization

DOCSTRING

# Fields:
- `orientation_angle::Float64`: DESCRIPTION
"""
struct LinearPolarization <: AbstractPolarization
    orientation_angle::Float64
end

"""
    ScalarPolarization

"Fake" polarization used to represent a scalar wave.
"""
struct ScalarPolarization <: AbstractPolarization end
struct RandomPolarization <: AbstractPolarization end
struct CircularPolarization <: AbstractPolarization end
struct EllipticalPolarization <: AbstractPolarization
    orientation_angle::Float64
    ellipticity_angle::Float64
end

abstract type AbstractField end


"""
    complex_amplitude(f::AbstractField)

TBW
"""
complex_amplitude(f::AbstractField) = error("Amplitude for $(typeof(f)) is not defined")
polarisation(f::AbstractField) = error("Polarisation for $(typeof(f)) is not defined")
wavelength(f::AbstractField) = error("Wavelength for $(typeof(f)) is not defined")
domain(f::AbstractField) = error("Domain for $(typeof(f)) is not defined")
# collect(f::AbstractField) = error("Collect for $(typeof(f)) is not defined")


struct SampledField{TA,TP} <: AbstractField
    amplitude::TA
    polarization::TP
    wavelength::Float64
    domain::CartesianDomain2D
end

complex_amplitude(f::SampledField) = f.amplitude
polarisation(f::SampledField) = f.polarization
wavelength(f::SampledField) = f.wavelength
domain(f::SampledField) = f.domain
# collect(f::SampledField) =
# field(amplitude(complex_amplitude(f)), phase(complex_amplitude(f)))

# TODO add multiplication of a field by a complex array
# import Base.*


struct LaserBeam{TP}
    wavelength::Float64
    polarization::TP
end
