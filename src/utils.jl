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
    rescale(array)  

Rescale array between 0 and 1
"""
function rescale(array)
    amin, amax = extrema(array)
    s = amin ≈ amax ? 1 : amax-amin
    return (array .- amin) / s
end



"""
    logrescale(array, α = 5)

Rescale in log-scalse: maximum will correspond to 1, ``10^{-α}`` will correspond to 0.
"""
function logrescale(array, α = 5)
    amin, amax = extrema(array)
    threshold = amax * 10. ^(-α)
    small = array .< threshold
    ret = copy(array)
    ret[small] .= threshold
    return rescale(log.(10,ret))
end

"Wrap Phase" phwrap(x::Real) = isnan(x) ? NaN : rem2pi(x, RoundNearest)
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
    ap = [ ((xc-o[1])^2 + (yc-o[2])^2) <= d^2 /4 ? 1 : 0 for yc ∈ yrange,  xc ∈ xrange]
    # area = +(ap[:]...)
    phmask = [ ((xc-o[1])^2 + (yc-o[2])^2) <= d^2 /4 ? 1 : NaN for yc ∈ yrange,  xc ∈ xrange]
    return(ap, phmask)
end

aperture(dom::CartesianDomain2D, d, o=(0,0)) = aperture(dom.xrange, dom.yrange, d,o)


"""
subdivide_sum(arr,Q)

Divide `arr`ay in `Q × Q` cells and sum the elements with the same indexes in each cell.
"""
function subdivide_sum(arr,Q)
    m,n = size(arr) .÷ Q
    ret = zeros(eltype(arr), m,n)

        for q in 0:Q^2-1
            a = q%Q
            b = q÷Q
            ret += @view arr[(m*a + 1) : m*(a+1) , (n*b + 1) : n*(b+1)]
        end

    return ret
end 