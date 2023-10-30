"""
Big class of general phase-retrieval-related problems.
"""
abstract type AbstractPRproblem <: Problem end

"""
  PRproblem(a, A)

Classical phase-retrieval problem of finding complex arrays x, X such that
  |x| = a, |X| = A, and X =F(x),
where F denotes the Fourier transform.
"""
struct PRproblem{T<:Real,N} <: AbstractPRproblem
    a::Array{T,N}
    A::Array{T,N}
end

# this problem can be converted to a feasibility problem
function TwoSetsFP(pr::PRproblem)
    return TwoSetsFP(
        ConstrainedByAmplitude(pr.a), FourierTransformedSet(ConstrainedByShape(pr.A))
    )
end

# methods to solve
"""
Gerchberg-Saxton, classical method to solve PR problem
"""
abstract type GS <: IterativeAlgorithm end

"""
    GSparam(; ϕ⁰, maxit, maxϵ, keephistory, snapshots)

Parameters for a concrete realisation of GS method.
"""
Base.@kwdef struct GSparam <: GS
    ϕ⁰ = missing
    maxϵ::Union{Float64,Missing} = missing
    maxit::Union{Missing,Int64} = missing
    keephistory::Bool = false
    snapshots::Array{Int64} = Int64[]
end

initial(alg::GSparam) = alg.ϕ⁰
tolerance(alg::GSparam) = alg.maxϵ
maxit(alg::GSparam) = alg.maxit
keephistory(alg::GSparam) = alg.keephistory
snapshots(alg::GSparam) = alg.snapshots

# GS is just alternating projection method
"""
    solve(pr::PRproblem, alg::GS, args...)

Solve PR problem using GS method. Any change of the GS parameters can be given
"""
function solve(pr::PRproblem, ::GSparam, ϕ⁰, maxϵ, maxit, keephistory, snapshots)
    return solve(
        TwoSetsFP(pr),
        AlternatingProjections.APparam(),
        ismissing(ϕ⁰) ? missing : pr.a .* exp.(2π .* im .* ϕ⁰),
        maxϵ,
        maxit,
        keephistory,
        snapshots,
    )
end

# PRproblem(x, X) = TwoSetsFP(ConstrainedByAmplitude(x), FourierTransformedSet(ConstrainedByAmplitude(X)))

"""
appsftoPR(ap,psfimage) constructs a phase retrieval problem from two real arrays, representing the pupil and the focal planes intensity
distributions.

The aperture and PSF are assumed to be centred in the array.

Examples
====
```jldoctest

ap,_ = PhaseRetrieval.aperture(-1:.2:1,-1:.2:1,.8);
psfimage = psf(ap);

pr = appsftoPR(ap,psfimage)

# output
PRproblem{Float64, 2}([1.0 0.0 … 1.0 1.0; 0.0 0.0 … 0.0 1.0; … ; 1.0 0.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0], [8.857504209336042 4.432323038895129 … 8.857504209336042 10.878351222990858; 4.432323038895129 0.7324956483358351 … 4.432323038895129 6.182768610120748; … ; 8.857504209336042 4.432323038895129 … 8.857504209336042 10.878351222990858; 10.878351222990858 6.182768610120747 … 10.878351222990858 12.999999999999998])
```

"""
function appsftoPR(ap, psfimage)
    size(ap) == size(psfimage) || error("Array sizes do not match")
    a = sqrt.(fftshift(ap))

    A = sqrt.(fftshift(psfimage))
    # normalise A to ap independent from FFT definition
    A = A ./ sqrt(sum(abs2, A)) .* sqrt(sum(abs2, fft(a)))

    return PRproblem(a, A)
end
