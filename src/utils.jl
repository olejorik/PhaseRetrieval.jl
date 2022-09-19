
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
    return kx*x .+ ky*y .+ k0
end

function lineararray(xrange::AbstractRange, yrange::AbstractRange, kx, ky, k0=0)
    return  kx*xrange' .+ ky*yrange .+ k0
end

lineararray(xrange::AbstractRange, yrange::AbstractRange, a::Vector, k0=0) = lineararray(xrange,yrange,a...,k0)

"""
    slopearray(dom::CartesianDomain2D, dx, dy, k0=0) generates linear array defined on the domain such that
    the total growth in x and y direction is dx and dy.
"""
function slopearray(dom::CartesianDomain2D, dx, dy, k0=0)
    xrange = dom.xrange
    yrange = dom.yrange
    kx = dx/ xrange.len
    ky = dy/ yrange.len
    return lineararray(xrange, yrange, kx, ky, k0)
end

function centroid(array)
    m00 = sum(array)
    m10 = sum([i*v for (i,v) in enumerate(sum(array, dims=1))])
    m01 = sum([i*v for (i,v) in enumerate(sum(array, dims=2))])
    return (m10/m00, m01/m00)
end

"""
    rescale(array)  

Rescale array between 0 and 1
"""
function rescale(array)
    amin, amax = extrema(array)
    s = amin ≈ amax ? 1 : amax-amin
    return (array .- amin) / Float64(s)
end

"""
    rescale_minmax(array)  

Rescale array between 0 and 1 and return the rescaled array and the min, max values
"""
function rescale_minmax(array)
    amin, amax = extrema(array)
    s = amin ≈ amax ? 1 : amax-amin
    return (array .- amin) / s, amin, amax
end



"""
    logrescale(array, α = 5)

Rescale in log-scale: maximum will correspond to 1, ``10^{-α}`` will correspond to 0.
"""
function logrescale(array, α = 5)
    amin, amax = extrema(array)
    threshold = amax * 10. ^(-α)
    small = array .< threshold
    ret = copy(array)
    ret[small] .= threshold
    return rescale(log10.(ret))
end

"Wrap Phase" phwrap(x::AbstractFloat) = isnan(x) ? NaN : rem2pi(x, RoundNearest)
# phwrap(::Val{NaN}) = NaN

function rotationmatrix(α)
    [cos(α) -sin(α); sin(α) cos(α)]
end


"""
    Generate quadratic array
"""
function quadratic(arrsize::Tuple)
    middle = (1 .+ arrsize)./2
    return reshape(1-middle[1]:arrsize[1]-middle[1], (:,1)).^2 .+ reshape(1-middle[2]:arrsize[2]-middle[2], (1,:)).^2;
end

quadratic(arrsize::Real) = quadratic((arrsize, arrsize))

"""
    Create array of zeroes with disk of ones of the relative diameter r
"""
function diskmatrix(gridsize::Integer, r=1)
    x = range(-1, 1, length=gridsize)
    y = range(-1, 1, length=gridsize)
    δ = 0. # tuning of the aperture size
    rd =r + δ  /gridsize 
    ap = [ (xc^2 + yc^2) <= r^2 ? 1 : 0 for  yc ∈ y, xc ∈ x]
    # area = +(ap[:]...)
    phmask = [ (xc^2 + yc^2) <= r^2 ? 1 : NaN for  yc ∈ y, xc ∈ x]
    return(ap, phmask)
end


"""
    aperture(xrange, yrange, d, o=0)
    aperture(dom::CartesianDomain2D, d, o)

Create circular aperture in array `xrange × yrange` with diameter `d` and center at `o`.
"""


function aperture(xrange::AbstractRange, yrange::AbstractRange, d, o=(0,0))
    ap = [ ((xc-o[1])^2 + (yc-o[2])^2) <= d^2 /4 ? 1. : 0. for yc ∈ yrange,  xc ∈ xrange]
    # area = +(ap[:]...)
    phmask = [ ((xc-o[1])^2 + (yc-o[2])^2) <= d^2 /4 ? 1. : NaN for yc ∈ yrange,  xc ∈ xrange]
    return(ap, phmask)
end

aperture(dom::CartesianDomain2D, d, o=(0,0)) = aperture(dom.xrange, dom.yrange, d,o)


"""
subdivide_sum(arr,Q)

Divide `arr`ay in quadratic cells of size `Q × Q`  and sum the elements with the same indexes in each cell.
"""
function subdivide_sum(arr,Q)
    if Q == 1
        return arr
    else
        m,n = size(arr) .÷ Q
        ret = zeros(eltype(arr), m,n)

            for q in 0:Q^2-1
                a = q%Q
                b = q÷Q
                ret += @view arr[(m*a + 1) : m*(a+1) , (n*b + 1) : n*(b+1)]
            end

        return ret
    end
end 

"""
    tile(arr, Q)

Divide `arr`ay in quadratic cells of size `Q × Q` and stack them along the third dimension.
"""
function tile(arr, Q :: Integer)    
    m,n = size(arr) .÷ Q
    B = reshape(arr, (Q,m,Q,n))
    C = PermutedDimsArray(B,[1,3,2,4])
    D = reshape(C, (Q,Q,:))
    
end

function tile(arr, Q :: Tuple)    
    m,n = size(arr) .÷ Q # TODO rewrite for any dimension
    B = reshape(arr, (Q[1],m,Q[2],n))
    C = PermutedDimsArray(B,[1,3,2,4])
    D = reshape(C, (Q[1],Q[2],:))
    
end

"""
    binning(arr, Q)

Downsample  `arr`ay by replacing quadratic cells of size `Q × Q` by summ of its elements
"""
function binning(arr, Q :: Integer)    
    m,n = size(arr) .÷ Q
    B = reshape(arr, (Q,m,Q,n))
    reshape(sum(B, dims=(1,3)),m,n) ./ Q^2
end


"""
    subdivide(arr, Q)

Subdivide array of dimension MQ x NQ in Q^2 stacked tiles of size M x N.
"""
function subdivide(arr, Q)
    Q1 = size(arr) .÷ Q 
    return tile(arr, Q1) 
end

function upscaleFactor(s, f, apertureD, λ )
    q = apertureD* s /(f *λ)
    upscale = ceil(Int,q)
end

function upscaleFactor(ims::ImagingSensor, λ )
    upscaleFactor(ims.cam.pixelsize, ims.lens.focallength, ims.lens.aperture, λ )
end

make_centered_domain2D(ims::ImagingSensor) = make_centered_domain2D(ims.cam.imagesize..., ims.cam.pixelsize)




"""
appsftoPR(ap,psfimage) constructs a phase retrieval problem from two real arrays, representing the pupil and the focal planes intensity
distributions. 

The aperture and PSF are assumed to be centred in the array.

Examples
====
```jldoctest

julia> ap,_ = PhaseRetrieval.aperture(-1:.2:1,-1:.2:1,.8)
([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 1.0 1.0 1.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0], [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN 1.0 NaN NaN NaN NaN NaN; NaN NaN NaN NaN 1.0 1.0 1.0 NaN NaN NaN NaN; NaN NaN NaN 1.0 1.0 1.0 1.0 1.0 NaN NaN NaN; NaN NaN NaN NaN 1.0 1.0 1.0 NaN NaN NaN NaN; NaN NaN NaN NaN NaN 1.0 NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN])

julia> psfimage = psf(ap)
[17.720305330423212 17.720305330423212 2.830830026003772 0.7990467619253477 1.7153703234534299 0.002318499057199975 0.8566413660014238 0.0023184990571999695 1.7153703234534299 0.7990467619253475 2.830830026003772; 17.720305330423212 17.720305330423212 2.830830026003772 0.7990467619253475 1.7153703234534299 0.0023184990571999695 0.8566413660014238 0.002318499057199975 1.7153703234534299 0.7990467619253477 2.830830026003772; 2.830830026003773 2.830830026003773 0.2240431494891658 5.881503709461158 4.671643508435733 0.0810140527710051 0.6181197482998197 0.0810140527710051 4.671643508435733 5.881503709461155 0.2240431494891658; 0.7990467619253475 0.7990467619253477 5.881503709461158 11.063720826850961 3.6825070656623606 0.6902785321094302 4.96008586865757 0.6902785321094302 3.6825070656623606 11.063720826850961 5.881503709461158; 1.7153703234534294 1.7153703234534299 4.671643508435735 3.6825070656623615 0.5365498748309357 19.64548752112056 38.22662768629444 19.64548752112056 0.5365498748309356 3.682507065662362 4.671643508435733; 0.00231849905719993 0.002318499057199919 0.08101405277100532 0.6902785321094301 19.64548752112056 78.45538081840574 118.33852533074673 78.45538081840574 19.64548752112056 0.6902785321094301 0.08101405277100517; 0.8566413660014239 0.8566413660014239 0.6181197482998197 4.96008586865757 38.22662768629444 118.33852533074673 169.0 118.33852533074673 38.22662768629444 4.96008586865757 0.6181197482998197; 0.002318499057199919 0.00231849905719993 0.08101405277100517 0.6902785321094301 19.64548752112056 78.45538081840574 118.33852533074673 78.45538081840574 19.64548752112056 0.6902785321094301 0.08101405277100532; 1.7153703234534299 1.7153703234534294 4.671643508435733 3.682507065662362 0.5365498748309356 19.64548752112056 38.22662768629444 19.64548752112056 0.5365498748309357 3.6825070656623615 4.671643508435735; 0.7990467619253477 0.7990467619253475 5.881503709461158 11.063720826850961 3.6825070656623606 0.6902785321094302 4.96008586865757 0.6902785321094302 3.6825070656623606 11.063720826850961 5.881503709461158; 2.830830026003773 2.830830026003773 0.2240431494891658 5.881503709461155 4.671643508435733 0.0810140527710051 0.6181197482998197 0.0810140527710051 4.671643508435733 5.881503709461158 0.2240431494891658]

julia> pr = appsftoPR(ap,psfimage);
AlternatingProjections.TwoSetsFP(AlternatingProjections.ConstrainedByAmplitude{Float64, 2}([1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0; 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0]), AlternatingProjections.FourierTransformedSet{AlternatingProjections.ConstrainedByAmplitude{Float64, 2}, FFTW.cFFTWPlan{ComplexF64, -1, false, 2, UnitRange{Int64}}, AbstractFFTs.ScaledPlan{ComplexF64, FFTW.cFFTWPlan{ComplexF64, 1, true, 2, UnitRange{Int64}}, Float64}}(AlternatingProjections.ConstrainedByAmplitude{Float64, 2}([13.0 10.878351222990862 6.182768610120748 2.2271250231312947 0.7862059197817196 0.9255492239753778 0.9255492239753778 0.7862059197817196 2.2271250231312947 6.182768610120748 10.878351222990862; 10.878351222990862 8.857504209336044 4.43232303889513 0.830830026003773 0.28462967654657045 0.048150794979936934 0.04815079497993703 0.2846296765465702 0.830830026003773 4.43232303889513 8.857504209336044; 6.182768610120748 4.43232303889513 0.7324956483358354 1.9189859472289947 2.1613985075491597 1.3097214678905702 1.30972146789057 2.1613985075491593 1.9189859472289947 0.7324956483358352 4.43232303889513; 2.2271250231312947 0.8308300260037731 1.9189859472289943 3.3262171947801247 2.425181170440913 0.893894155884995 0.8938941558849948 2.425181170440913 3.3262171947801247 1.9189859472289943 0.8308300260037731; 0.7862059197817196 0.28462967654657007 2.1613985075491593 2.425181170440913 0.4733319654208511 1.6825070656623624 1.6825070656623624 0.4733319654208511 2.4251811704409127 2.1613985075491593 0.28462967654657007; 0.9255492239753776 0.04815079497993745 1.3097214678905702 0.8938941558849948 1.6825070656623622 4.209549302529098 4.209549302529098 1.6825070656623622 0.893894155884995 1.3097214678905702 0.04815079497993751; 0.9255492239753776 0.04815079497993751 1.3097214678905702 0.893894155884995 1.6825070656623622 4.209549302529098 4.209549302529098 1.6825070656623622 0.8938941558849948 1.3097214678905702 0.04815079497993745; 0.7862059197817196 0.28462967654657007 2.1613985075491593 2.4251811704409127 0.4733319654208511 1.6825070656623624 1.6825070656623624 0.4733319654208511 2.425181170440913 2.1613985075491593 0.28462967654657007; 2.2271250231312947 0.8308300260037731 1.9189859472289943 3.3262171947801247 2.425181170440913 0.8938941558849948 0.893894155884995 2.425181170440913 3.3262171947801247 1.9189859472289943 0.8308300260037731; 6.182768610120748 4.43232303889513 0.7324956483358352 1.9189859472289947 2.1613985075491593 1.30972146789057 1.3097214678905702 2.1613985075491597 1.9189859472289947 0.7324956483358354 4.43232303889513; 10.878351222990862 8.857504209336044 4.43232303889513 0.830830026003773 0.2846296765465702 0.04815079497993703 0.048150794979936934 0.28462967654657045 0.830830026003773 4.43232303889513 8.857504209336044]), FFTW forward plan for 11×11 array of ComplexF64
(dft-rank>=2/1
  (dft-direct-11-x11 "n1fv_11_avx2_128")
  (dft-direct-11-x11 "n1fv_11_avx2_128")), 0.008264462809917356 * FFTW in-place backward plan for 11×11 array of ComplexF64
(dft-rank>=2/1
  (dft-direct-11-x11 "n1bv_11_avx2_128")
  (dft-direct-11-x11 "n1bv_11_avx2_128"))))

```

"""
function appsftoPR(ap,psfimage)
    size(ap) == size(psfimage) || error("Array sizes do not match")
    a = sqrt.(fftshift(ap))

    A = sqrt.(fftshift(psfimage))
    # normalise A to ap independent from FFT definition
    A = A ./ sqrt(sum(abs2,A)) .* sqrt(sum(abs2,fft(a)))

    return PRproblem(a, A)

end

function removepiston(ϕ)
    return  ϕ  .- mean(phwrap.(filter(!isnan,ϕ)))
end

function removepiston(ϕ, mask)
    return  ϕ  .- mean(phwrap.(filter(!isnan, ϕ .* mask)))
end

function removetiptilt(ϕ)
    sy, sx = size(ϕ)
    dx = diff(ϕ, dims=2)
    dy = diff(ϕ, dims=1)
    kx = mean(phwrap.(filter(!isnan,dx)))
    ky = mean(phwrap.(filter(!isnan,dy)))
    tiptilt = lineararray(1:sx,1:sy, kx,ky)
    return  phwrap.(ϕ  .- tiptilt)
end

function removetiptilt(ϕ, mask)
    sy, sx = size(ϕ)
    dx = diff(ϕ .* mask, dims=2)
    dy = diff(ϕ .* mask, dims=1)
    kx = mean(phwrap.(filter(!isnan,dx)))
    ky = mean(phwrap.(filter(!isnan,dy)))
    tiptilt = lineararray(1:sx,1:sy, kx,ky)
    return  phwrap.(ϕ  .- tiptilt)
end

function twinphase(ϕ)
    shifts = map(x->1-mod(x,2), size(ϕ))
    return circshift(reverse(-ϕ),shifts)
end