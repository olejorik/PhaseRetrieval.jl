# Different methods of generating PSF

using FFTW

field(amplitude, phase) = amplitude .* exp.( 1im * phase)

"Calculate psf for given amplitude and phase"
psf(field::Array) = abs.(toimageplane(field)) .^2

psf(amplitude, phase) = psf(field(amplitude, phase))

subpsf(amplitude, phase,Q) = psf(subdivide_sum(field(amplitude, phase),Q))



toimageplane(field) = ifftshift(fft(fftshift(field)))



