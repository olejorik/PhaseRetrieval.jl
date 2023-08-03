
"""
    lineararray(xrange, yrange, kx, ky, k0=0)

Evaluate values of linear function ``kx⋅x + ky⋅y + k0``  on array  ``xrange × yrange``.

    lineararray(size, kx, ky, k0=0)

Use `1:size` as `xrange` and `yrange`.

    lineararray(xrange, yrange, a::Vector, k0=0)

Use first and second components of vector `a` as `kx` and `ky`.

"""
function lineararray(size, kx, ky, k0=0)
    x = (1:size)'
    y = 1:size
    return kx * x .+ ky * y .+ k0
end

function lineararray(xrange::AbstractRange, yrange::AbstractRange, kx, ky, k0=0)
    return kx * xrange' .+ ky * yrange .+ k0
end

function lineararray(xrange::AbstractRange, yrange::AbstractRange, a::Vector, k0=0)
    return lineararray(xrange, yrange, a..., k0)
end

"""
    slopearray(dom::CartesianDomain2D, dx, dy, k0=0) generates linear array defined on the domain such that
    the total growth in x and y direction is dx and dy.
"""
function slopearray(dom::CartesianDomain2D, dx, dy, k0=0)
    xrange = dom.xrange
    yrange = dom.yrange
    kx = dx / xrange.len
    ky = dy / yrange.len
    return lineararray(xrange, yrange, kx, ky, k0)
end

function centroid(array)
    m00 = sum(array)
    m10 = sum([i * v for (i, v) in enumerate(sum(array; dims=1))])
    m01 = sum([i * v for (i, v) in enumerate(sum(array; dims=2))])
    return (m10 / m00, m01 / m00)
end

"""
    rescale(array)

Rescale array between 0 and 1
"""
function rescale(array)
    amin, amax = extrema(array)
    s = amin ≈ amax ? 1 : amax - amin
    return (array .- amin) / Float64(s)
end

"""
    rescale_minmax(array)

Rescale array between 0 and 1 and return the rescaled array and the min, max values
"""
function rescale_minmax(array)
    amin, amax = extrema(array)
    s = amin ≈ amax ? 1 : amax - amin
    return (array .- amin) / s, amin, amax
end

"""
    logrescale(array, α = 5)

Rescale in log-scale: maximum will correspond to 1, ``10^{-α}`` will correspond to 0.
"""
function logrescale(array, α=5)
    amin, amax = extrema(array)
    threshold = amax * 10.0^(-α)
    small = array .< threshold
    ret = copy(array)
    ret[small] .= threshold
    return rescale(log10.(ret))
end

"Wrap Phase"
phwrap(x::AbstractFloat) = isnan(x) ? NaN : rem2pi(x, RoundNearest)
# phwrap(::Val{NaN}) = NaN

function rotationmatrix(α)
    return [cos(α) -sin(α); sin(α) cos(α)]
end

"""
    Generate quadratic array
"""
function quadratic(arrsize::Tuple)
    middle = (1 .+ arrsize) ./ 2
    return reshape((1 - middle[1]):(arrsize[1] - middle[1]), (:, 1)) .^ 2 .+
           reshape((1 - middle[2]):(arrsize[2] - middle[2]), (1, :)) .^ 2
end

quadratic(arrsize::Real) = quadratic((arrsize, arrsize))

"""
    Create array of zeroes with disk of ones of the relative diameter r
"""
function diskmatrix(gridsize::Integer, r=1)
    x = range(-1, 1; length=gridsize)
    y = range(-1, 1; length=gridsize)
    δ = 0.0 # tuning of the aperture size
    rd = r + δ / gridsize
    ap = [(xc^2 + yc^2) <= r^2 ? 1 : 0 for yc in y, xc in x]
    # area = +(ap[:]...)
    phmask = [(xc^2 + yc^2) <= r^2 ? 1 : NaN for yc in y, xc in x]
    return (ap, phmask)
end

"""
    aperture(xrange, yrange, d, o=0)
    aperture(dom::CartesianDomain2D, d, o)

Create circular aperture in array `xrange × yrange` with diameter `d` and center at `o`.
"""

function aperture(xrange::AbstractRange, yrange::AbstractRange, d, o=(0, 0))
    ap = [
        ((xc - o[1])^2 + (yc - o[2])^2) <= d^2 / 4 ? 1.0 : 0.0 for yc in yrange,
        xc in xrange
    ]
    # area = +(ap[:]...)
    phmask = [
        ((xc - o[1])^2 + (yc - o[2])^2) <= d^2 / 4 ? 1.0 : NaN for yc in yrange,
        xc in xrange
    ]
    return (ap, phmask)
end

aperture(dom::CartesianDomain2D, d, o=(0, 0)) = aperture(dom.xrange, dom.yrange, d, o)

"""
subdivide_sum(arr,Q)

Divide `arr`ay in quadratic cells of size `Q × Q`  and sum the elements with the same indexes in each cell.
"""
function subdivide_sum(arr, Q)
    if Q == 1
        return arr
    else
        m, n = size(arr) .÷ Q
        ret = zeros(eltype(arr), m, n)

        for q in 0:(Q^2 - 1)
            a = q % Q
            b = q ÷ Q
            ret += @view arr[(m * a + 1):(m * (a + 1)), (n * b + 1):(n * (b + 1))]
        end

        return ret
    end
end

"""
    tile(arr, Q)

Divide `arr`ay in quadratic cells of size `Q × Q` and stack them along the third dimension.
"""
function tile(arr, Q::Integer)
    m, n = size(arr) .÷ Q
    B = reshape(arr, (Q, m, Q, n))
    C = PermutedDimsArray(B, [1, 3, 2, 4])
    return D = reshape(C, (Q, Q, :))
end

function tile(arr, Q::Tuple)
    m, n = size(arr) .÷ Q # TODO rewrite for any dimension
    B = reshape(arr, (Q[1], m, Q[2], n))
    C = PermutedDimsArray(B, [1, 3, 2, 4])
    return D = reshape(C, (Q[1], Q[2], :))
end

"""
    binning(arr, Q)

Downsample  `arr`ay by replacing quadratic cells of size `Q × Q` by summ of its elements
"""
function binning(arr, Q::Integer)
    m, n = size(arr) .÷ Q
    B = reshape(arr, (Q, m, Q, n))
    return reshape(sum(B; dims=(1, 3)), m, n) ./ Q^2
end

"""
    subdivide(arr, Q)

Subdivide array of dimension MQ x NQ in Q^2 stacked tiles of size M x N.
"""
function subdivide(arr, Q)
    Q1 = size(arr) .÷ Q
    return tile(arr, Q1)
end

function ap_ratio(s, f, apertureD, λ)
    q = apertureD * s / (f * λ)
    return q
end

function ap_ratio(ims::ImagingSensor, λ)
    return ap_ratio(ims.cam.pixelsize, ims.lens.focallength, ims.lens.aperture, λ)
end

function ap_ratio(c::SimConfig)
    return ap_ratio(c.ims, c.λ)
end

function upscaleFactor(s, f, apertureD, λ)
    return ceil(Int, ap_ratio(s, f, apertureD, λ))
end

function upscaleFactor(ims::ImagingSensor, λ)
    return upscaleFactor(ims.cam.pixelsize, ims.lens.focallength, ims.lens.aperture, λ)
end

function upscaleFactor(c::SimConfig)
    return upscaleFactor(c.ims, c.λ)
end

function make_centered_domain2D(ims::ImagingSensor)
    return make_centered_domain2D(ims.cam.imagesize..., ims.cam.pixelsize)
end
function make_centered_domain2D(hw::hwConfig)
    return make_centered_domain2D(hw.cam.imagesize..., hw.cam.pixelsize)
end

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

function removepiston(ϕ)
    return ϕ .- mean(phwrap.(filter(!isnan, ϕ)))
end

function removepiston(ϕ, mask)
    return ϕ .- mean(phwrap.(filter(!isnan, ϕ .* mask)))
end

function removetiptilt(ϕ)
    sy, sx = size(ϕ)
    dx = diff(ϕ; dims=2)
    dy = diff(ϕ; dims=1)
    kx = mean(phwrap.(filter(!isnan, dx)))
    ky = mean(phwrap.(filter(!isnan, dy)))
    tiptilt = lineararray(1:sx, 1:sy, kx, ky)
    return phwrap.(ϕ .- tiptilt)
end

function removetiptilt(ϕ, mask)
    sy, sx = size(ϕ)
    dx = diff(ϕ .* mask; dims=2)
    dy = diff(ϕ .* mask; dims=1)
    kx = mean(phwrap.(filter(!isnan, dx)))
    ky = mean(phwrap.(filter(!isnan, dy)))
    tiptilt = lineararray(1:sx, 1:sy, kx, ky)
    return phwrap.(ϕ .- tiptilt)
end

function twinphase(ϕ)
    shifts = map(x -> 1 - mod(x, 2), size(ϕ))
    return circshift(reverse(-ϕ), shifts)
end
