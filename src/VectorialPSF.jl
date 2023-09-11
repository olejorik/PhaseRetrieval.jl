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
