# The Field structure should represent the optical field in a plane. It should contain the values of the complex amplitude, wavelength, and the polarisation.
#

export ComplexAmplitude,
    ScalarPolarization,
    RandomPolarization,
    CircularPolarization,
    EllipticalPolarization,
    Field,
    LaserBeam


abstract type AbstractComplexAmplitude end
struct ComplexAmplitude{TA,TPH} <: AbstractComplexAmplitude
    amplitude::TA
    phase::TPH
end

abstract type AbstractPolarization end
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
    elliptcity_angle::Float64
end

abstract type AbstractField end

struct Field{TA,TP} <: AbstractField
    amplitude::TA
    polarization::TP
    wavelength::Float64
end

struct LaserBeam{TP}
    wavelength::Float64
    polarization::TP
end
