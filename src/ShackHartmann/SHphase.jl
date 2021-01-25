using FFTW


# Algorithms to choose
abstract type Alg end
"Build SH diversity using the quadratic functions" struct Quadratic <: Alg end
"Build SH diversity using the linear functions" struct Linear <: Alg end
   
using ..PhaseRetrieval: lineararray

"""
    SHdiversity(xrange, yrange, e₁, e₂, o)

Generate SH diversity phase on array `xrange × yrange` with cell based on orthogonal grid
formed by vectors  `e₁, e₂,` and origin at `o`. 
"""
function SHdiversity(size, cellsize, celloffset=0, α=0)
    shd = zeros(size)
    
    # # first generate indexes
    # ind = orthIndexes(size, cellsize, celloffset, α)
    
    # # Now find range of indexes
    # irange = extrema(ind[:,:,1])
    # jrange = extrema(ind[:,:,2])
    # indoffset = -first.([irange,jrange]) .+ 1
    
    # # generate array of cell centres
    # centres = [centre(i, j, cellsize, celloffset, α) 
    #             for i in UnitRange(irange...), j in UnitRange(jrange...)]
     
    (ind, centres, indoffset) = indexes_centres(size, cellsize, celloffset, α)

    # Method through MLA phase and global defocus
    # Now for every point in the array take the distance squared to centre of its cell
    # This will construct MLA phase
    for i in 1:size[1], j in 1:size[2]
        cellind = ind[i,j,:]
        cellcentre = centres[(cellind .+ indoffset)...]
        shd[i,j] = +((([i,j] - cellcentre).^2)...)
    end
        
        
    # Now we can subtract the global defocus
    defocus =  quadratic(size) 
    shd -= defocus
    
    return shd        
end

SHdiversity(size::Real, cellsize, celloffset, α) = SHdiversity((size, size), cellsize, celloffset, α)

function SHdiversity(size, cellsize, celloffset, α, ::Linear)
    shd = zeros(size)

    (ind, centres, indoffset) = indexes_centres(size, cellsize, celloffset, α)

    # Method through explicit piece-wise linear function
    # - 2 c_i x + c_i^2, with x = i -size/2
    for i in 1:size, j in 1:size
        cellind = ind[i,j,:]
        cellcentre = centres[(cellind .+ indoffset)...]
        shd[i,j] = -2*cellcentre .* (cellind .- size/2) .+ dot(cellcentre, cellcentre)
                
    end

    return shd        
end

SHdiversity(size::Real, cellsize, celloffset, α, method) = SHdiversity((size,size), cellsize, celloffset, α, method)

function SHdiversity(xrange::AbstractRange, yrange::AbstractRange, e₁, e₂, o)
    cellind = indicesmap(xrange, yrange, e₁, e₂,o)
    irange  = UnitRange(extrema(cellind[1])...)
    jrange  = UnitRange(extrema(cellind[2])...)

    cellcentres = centresarray(cellind, e₁, e₂,o)

    s0 = 0.5*(e₁ + e₂) .+ o # centre of the zeroth cell
    p1 = dot(e₁, e₁)
    p2 = dot(e₂, e₂)
    s01 = dot(e₁, s0)
    s02 = dot(e₂, s0)
    s00 = dot(s0, s0)

    ret = lineararray(xrange, yrange, 2*s0, s00) # l0 = (s0,2 x) +s0^2
    l1 = lineararray(xrange, yrange, e₁)
    l2 = lineararray(xrange, yrange, e₂)

    for i in irange
        zi = cellind[1] .== i
        ci = - i^2*p1  - 2i*s01 
        ret[zi] += 2i*l1[zi] .+ ci
    end

    for j in jrange
        zj = cellind[2] .== j
        cj = - j^2*p2  - 2j*s02 
        ret[zj] += 2j*l2[zj] .+ cj
    end
        



    return ret
end

function indicesmap(xrange, yrange, e₁, e₂,o)
    iind =fld.(lineararray(xrange, yrange, e₁, - dot(e₁, o)), dot(e₁, e₁))
    jind =fld.(lineararray(xrange, yrange, e₂, - dot(e₂, o)), dot(e₂, e₂))
    return(iind, jind)
    
end

function centresarray(cellind, e₁, e₂,o)
    ijrange=extrema.(cellind)
    cij = [ (i + 0.5)*e₁ + (j + 0.5)*e₂ .+  o  for i in UnitRange(ijrange[1]...), j in UnitRange(ijrange[2]...)  ]
    return cij   
end
 
# TODO finish SH diversity for a sensor and total size of the phase field
# function SHdiversity(wfs::SHSensor) = 
#     SHdiversity(wfs.CameraChip.imagesize, cellsize, celloffset=0, α=0)


"""
    Generate array of indexes of an orthogonal geometry
"""
function orthIndexes(arrsize::Tuple, cellsize, gridorigin=0, α=0)
    ind = zeros(arrsize...,2)
    R = rotationmatrix(-α)
    x = 1- gridorigin[1] : arrsize[1]- gridorigin[1]
    y = reshape(1 - gridorigin[2] : arrsize[2] - gridorigin[2], (1,:))
    # for i in 1:arrsize, j in 1:arrsize
    #     ind[:,i, j] = R * [i, j]
    # end
    ind[:,:,1] = x .* R[1,1] .+ y .* R[1,2]
    ind[:,:,2] = x .* R[2,1] .+ y .* R[2,2] 
    ind = fld.(ind , cellsize)
    return Int64.(ind)
end

orthIndexes(arrsize::Tuple, cellsize, gridorigin::Real, α=0) =
    orthIndexes(arrsize, cellsize, (gridorigin, gridorigin), α)

orthIndexes(arrsize::Real, cellsize, gridorigin, α=0) =
    orthIndexes((arrsize,arrsize), cellsize, gridorigin, α)

# this implementation appears to be quite slow
# function orthIndexes(arrsize, cellsize, gridorigin=0, α=0)
#     R = rotationmatrix(-α)
#     ind = [Int.(fld.(R * ( [i, j] .- gridorigin), cellsize)) 
#             for i in 1:arrsize, j in 1:arrsize]    
# end


"""
    Give centre of an orthogonal grid cell of the grid with cellsize, 
    centre at grid origin and rotated at angle α counterclockwise.
"""
function centre(i, j, cellsize, gridorigin, α)
    R = rotationmatrix(α)
    R * ([i,j] * cellsize  .+ cellsize / 2 ) .+ gridorigin
end

function indexes_centres(size, cellsize, celloffset, α)
    
    # first generate indexes
    ind = orthIndexes(size, cellsize, celloffset, α)
    
    # Now find range of indexes
    irange = extrema(ind[:,:,1])
    jrange = extrema(ind[:,:,2])
    indoffset = -first.([irange,jrange]) .+ 1
    
    # generate array of cell centres
    centres = [centre(i, j, cellsize, celloffset, α) 
                for i in UnitRange(irange...), j in UnitRange(jrange...)]
    
    return (ind, centres, indoffset)
end



"""
    Micro-lens array structure. 

    The main reason to introduce it separable is to handle the case when camera chip is located not exactly
    in the MLA focal plane.
"""
Base.@kwdef struct MLA
    pitch
    focallength
end

Base.@kwdef struct CameraChip
    pixelsize
    imagesize
end


"""
    Shack-Hartmann sensor composed from MLA and a camera. 
    
    MLA plane position is defined by distance,  β,  and γ. 
    MLA subapertures position is defined by mlaorigin and α.  
"""
Base.@kwdef struct SHSensor
    mla::MLA
    cam::CameraChip
    distance = mla.focallength
    mlaorigin = [0,0]
    α = 0
    β = 0
    γ = 0
end

struct HPPlan
    wfs::SHSensor
    width
    height

end

struct PFPlan
    xdomain::Domain{2}
    fdomain
    ap
    diversity
    # f
    # l
end




function upscaleFactor(wfs::SHSensor, apertureD, λ )
    s = wfs.cam.pixelsize
    f = wfs.mla.focallength
    q = apertureD* s /(f *λ)
    upscale = ceil(Int,q)
end

"""
    constructSHPhaseDiversity(wfs::SHSensor, d, λ )

Construct SH phase diversity of minimal size for `wfs` within aperture of diameter `d` and wavelength `λ`.
Return tuple of arrays for the aperture and phase.
"""
function constructSHPhaseDiversity(wfs::SHSensor, apertureD, λ )
    s = wfs.cam.pixelsize
    f = wfs.mla.focallength
    p = wfs.mla.pitch
    o = wfs.mlaorigin
    α = wfs.α   
    imsize = wfs.cam.imagesize

    # first check the dimensions
    q = apertureD* s /(f *λ)
    # if upscale=1, we are fine. otherwise we need to upscale
    upscale = upscaleFactor(wfs, apertureD, λ )
    cols,rows = imsize .* upscale
    # size of aperture
    arrsize = ceil(Int,apertureD/s * upscale)
    phasescale = pi/(1*arrsize) * q / upscale
    # or size of the camera chip
    # TODO extend to different scaling
    # arrsize =  imsize .* upscale
    # phasescale = pi/(1*arrsize[1]) * q / upscale
    cellsize = wfs.mla.pitch/s * upscale

    # in pixel coordianates
    # shd = SHdiversity(arrsize, cellsize, o, α) * phasescale
    # ap, apmask = diskmatrix(arrsize[1], q/upscale)

    # in physical coordinates
    # rng = -apertureD/2 : s/upscale : apertureD # wrong, HO size can be bigger
    # xrng= ((1:cols) .- ceil(cols/2)) * 1/(s*cols / upscale) * (λ * f)
    # yrng= ((1:rows) .- ceil(rows/2)) * 1/(s*rows / upscale) * (λ * f)

    # even better, through the domain
    wfsdom = make_centered_domain2D(imsize...,s)
    xrng = dualRange(wfsdom.xrange, upscale) * (λ * f)
    yrng = dualRange(wfsdom.yrange, upscale) * (λ * f)
    basis = rotationmatrix(α) * p 
    shd = - π / (f*λ) * SHdiversity(xrng, yrng, basis[:,1], basis[:,2],o  ) 
    ap, apmask = aperture(xrng,yrng, apertureD)
    return (ap,shd)
    
end

make_centered_domain2D(wfs::SHSensor) = make_centered_domain2D(wfs.cam.imagesize..., wfs.cam.pixelsize)
