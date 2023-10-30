# General functions for working with phase values

"Wrap Phase"
phwrap(x::AbstractFloat) = isnan(x) ? NaN : rem2pi(x, RoundNearest)
# phwrap(::Val{NaN}) = NaN

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
