# general functions for psf simulations

function ap_ratio(s, f, apertureD, λ)
    q = apertureD * s / (f * λ)
    return q
end

function upscaleFactor(s, f, apertureD, λ)
    return ceil(Int, ap_ratio(s, f, apertureD, λ))
end

"""
    aperture(xrange, yrange, d, o=0)
    aperture(dom::CartesianDomain2D, d, o)

Create circular aperture in array `xrange × yrange` with diameter `d` and center at `o`.
"""
function aperture(xrange::AbstractRange, yrange::AbstractRange, d, o=(0, 0))
    return diskmatrix(xrange, yrange, d; o=(0, 0))
end

aperture(dom::CartesianDomain2D, d, o=(0, 0)) = aperture(dom.xrange, dom.yrange, d, o)
