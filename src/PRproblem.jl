"""
Big class of general phase-retrieval-related problems.
"""
abstract type AbstractPRproblem <: Problem end

"""
  PRproblem(a, A)

Classical phase-retrieval problem of finding complex arrays x, X such that
  |x| = a, |X| = A, and X =F(x),
where F denotes the Fourier transform.
Note that the problem is not required to be consistent.
"""
struct PRproblem{T<:Real,N} <: AbstractPRproblem
    a::Array{T,N}
    A::Array{T,N}
end

"""
    prproblem(a,A; center = true, renormalize = true) -> PRproblem

Default constructor for PRproblem from the arrays which origin is supposed to be at their center and the second array is renormalised to satisfy the Parceval equation.
"""
function prproblem(a, A; center=true, renormalize=true)
    A1 = renormalize ? A ./ sqrt(sum(abs2, A)) .* sqrt(sum(abs2, a)) : A
    shift = center ? fftshift : x -> x
    return PRproblem(shift(a), shift(A1))
end

# PRproblem can be converted to a feasibility problem
function TwoSetsFP(pr::PRproblem)
    return TwoSetsFP(
        ConstrainedByAmplitude(pr.a), FourierTransformedSet(ConstrainedByShape(pr.A))
    )
end

# TODO add PRsat, add automatic ffshift and norming

"""
  PRproblemSat(a, A, satlevel)

Saturated PSF phase-retrieval problem of finding complex arrays x, X such that
  |x| = `a`, max(|X|, `satlevel`) = `A`, and X =F(x),
where F denotes the Fourier transform.
"""
struct PRproblemSat{T<:Real,N} <: AbstractPRproblem
    a::Array{T,N}
    A::Array{T,N}
    s::T
end

# this problem can be converted to a feasibility problem in two ways
function TwoSetsFP(pr::PRproblemSat, scaling="focal")
    sat = pr.A .< sqrt(pr.s)
    return TwoSetsFP(
        ConstrainedByAmplitude(pr.a),
        FourierTransformedSet(ConstrainedByShapeSaturated(pr.A, sat)),
    )
end

# methods to solve

abstract type PRalgorithm <: Algorithm end
"""
Gerchberg-Saxton, classical method to solve PR problem
"""
abstract type GS <: PRalgorithm end

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
#
APparam(alg::GSparam) =
    APparam(initial(alg), tolerance(alg), maxit(alg), keephistory(alg), snapshots(alg))
"""
    solve(pr::PRproblem, alg::GS, args...)

Solve PR problem using GS method. Any change of the GS parameters can be given
"""
solve(pr::PRproblem, alg::GS, args...; kwargs...) =
    solve(TwoSetsFP(pr), APparam(alg), args...; kwargs...)

solve(pr::PRproblem, alg::IterativeAlgorithm, args...; kwargs...) =
    solve(TwoSetsFP(pr), alg, args...; kwargs...)

solve(pr::PRproblem, algs::Tuple{Vararg{IterativeAlgorithm}}, args...; kwargs...) =
    solve(TwoSetsFP(pr), algs, args...; kwargs...)

# PRproblem(x, X) = TwoSetsFP(ConstrainedByAmplitude(x), FourierTransformedSet(ConstrainedByAmplitude(X)))
