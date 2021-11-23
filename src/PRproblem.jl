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
struct PRproblem{T <: Real, N} <: AbstractPRproblem
  a::Array{T,N}
  A::Array{T,N}
end

# this problem can be converted to a feasibility problem
TwoSetsFP(pr::PRproblem) = TwoSetsFP(ConstrainedByAmplitude(pr.a), FourierTransformedSet(ConstrainedByAmplitude(pr.A)))

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
    maxϵ::Union{Float64, Missing} = missing
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
solve(pr::PRproblem, ::GSparam, ϕ⁰, maxϵ, maxit, keephistory, snapshots) = 
        solve(
            TwoSetsFP(pr),            
            AlternatingProjections.APparam(),
            ismissing(ϕ⁰) ?  missing : pr.a .* exp.(2π .* im .* ϕ⁰),
            maxϵ,
            maxit,
            keephistory,
            snapshots)

# PRproblem(x, X) = TwoSetsFP(ConstrainedByAmplitude(x), FourierTransformedSet(ConstrainedByAmplitude(X)))
