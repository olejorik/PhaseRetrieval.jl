# Different methods of generating PSF

using FFTW

field(amplitude, phase) = amplitude .* exp.( 1im * phase)

"Calculate psf for given amplitude and phase"
psf(field::Array) = abs.(toimageplane(field)) .^2

psf(amplitude, phase) = psf(field(amplitude, phase))

subpsf(amplitude, phase,Q) = psf(subdivide_sum(field(amplitude, phase),Q))



toimageplane(field) = ifftshift(fft(fftshift(field)))

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
