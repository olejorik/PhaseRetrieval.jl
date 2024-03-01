export incoherent_psf

function get_cosines(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    σx, σy = getranges(ddom)
    σs = (σz(σ...) - 1 for σ in ddom)
    return (σx=σx, σy=σy, σz=σs)
end

function get_polarization_magnitudes(conf::SimConfig)
    ddom = conf.dualroi * (1 / focallength(conf))
    inaperture = (!isnan).(conf.mask)
    # σx, σy = getranges(ddom)
    exx = Array{Float64}(undef, size(ddom))
    eyx = Array{Float64}(undef, size(ddom))
    ezx = Array{Float64}(undef, size(ddom))
    exy = Array{Float64}(undef, size(ddom))
    eyy = Array{Float64}(undef, size(ddom))
    ezy = Array{Float64}(undef, size(ddom))
    for (i, (y, x)) in enumerate(ddom)
        if inaperture[i]
            zfactor = 1 / (1 + σz(x, y))
            exx[i] = 1 - x^2 * zfactor
            eyx[i] = -x * y * zfactor
            ezx[i] = -x
            exy[i] = -x * y * zfactor
            eyy[i] = 1 - x^2 * zfactor
            ezy[i] = -y
        end
    end
    return (exx=exx, eyx=eyx, ezx=ezx, exy=exy, eyy=eyy, ezy=ezy)
end

"""
    incoherent_psf(array<:Array)

Calculate 2D array corresponding to the incoherent sum of the PSF formed by slices `array[:,:,ind...]` for all possible indici `ind`.

```jldoctest
julia> incoherent_psf(ones(3,3,2,3))
3×3 Matrix{Float64}:
 0.0    0.0  0.0
 0.0  486.0  0.0
 0.0    0.0  0.0

julia> incoherent_psf(ones(4,4,2,2,1))
4×4 Matrix{Float64}:
 0.0  0.0     0.0  0.0
 0.0  0.0     0.0  0.0
 0.0  0.0  1024.0  0.0
 0.0  0.0     0.0  0.0

julia> incoherent_psf(ones(4,4))
4×4 Matrix{Float64}:
 0.0  0.0    0.0  0.0
 0.0  0.0    0.0  0.0
 0.0  0.0  256.0  0.0
 0.0  0.0    0.0  0.0

```
"""
function incoherent_psf(array::Array)
    arr_dims = ndims(array)
    arr_dims > 1 || error("Array should have at least 2 dimensions")

    ret = zeros(size(array)[1:2])
    for f in eachslice(array; dims=Tuple(3:arr_dims))
        ret .+= psf(f)
    end

    return ret
end
