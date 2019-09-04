#The "1/e^2 beamwidth" seems to have inconsistent definitions,
#either meaning the half-width or the full width of the beam
#Here "beamwidth" will refer to the full width, and "sigma" will
#refer to the halfwidth
#Sigma is the width of a gaussian defined as y = e^-(2x^2/sigma^2)
#(this parameterization of the Gaussian is common for optics but not universal. The other
#common convention has a 2 in the denominator of the exponent rather thatn the numerator.
#This is the convention in statistics and thus is how variance and std are defined.  So if we want
#an r.m.s. value we should use this latter convention)
#By the former convention "sigma" is the distance from the origin at which the gaussian
#decays to 1/e^2 of its maximum value, i.e. the halfwidth.
#sigma_statistics = std = sigma_optics/2
#sigma_optics = 2*std


#gaussian approximation to airy disk
#taken from paraxial approximations derived in:
#Gaussian approximations of fluorescence microscope point-spread function models
#Bo Zhang, Josiane Zerubia, and Jean-Christophe Olivo-Marin
#1 April 2007 Vol. 46, No. 10 APPLIED OPTICS
#Note that the authors derive the below expression for "sigma"
#but their definition of sigma differs from ours; they use the statistical convention
#(see comments at top of this file) instead of 1/e^2 convention
#their r.h.s. is    e^-(x^2/2*sigma^2)
#...while ours is   e^-(2x^2/sigma^2)
#thus our axial sigma is twice theirs
approx_sigma_axial(lambda, NA, n) = 4*sqrt(6)*n/(2*pi/lambda*NA^2)
approx_sigma_axial(lambda, f, unfocused_sigma, n) = approx_sigma_axial(lambda, lens_na(unfocused_sigma, f), n)

approx_sigma_lateral(lambda, NA) = sqrt(2)*lambda / (pi*NA)
approx_sigma_lateral(lambda, f, unfocused_sigma) = approx_sigma_lateral(lambda, lens_na(unfocused_sigma, f))
focused_beamwidth(lambda, f, unfocused_beamwidth) = 2*approx_sigma_lateral(lambda, f, unfocused_beamwidth/2)

#assuming a lens has a big enough aperture to accept
#a collimated input of specified halfwidth and also
#that the medium is the same before and after focusing
#what is the NA of the system?
lens_na(unfocused_sigma, f) = unfocused_sigma/f

#Halfwidth of the Gaussian approximation of the axial PSF of a light sheet imaging system,
#considering both the PSF of the imaging lens and the sheet thickness
lightsheet_axial_combined_psf(det_sigma, sheet_sigma) =
    sqrt((det_sigma^2*sheet_sigma^2)/(det_sigma^2 + sheet_sigma^2))

rayleigh_length(lambda, focused_sigma) = pi*focused_sigma^2 / lambda
rayleigh_length(lambda, f, unfocused_sigma, n) =
    rayleigh_length(lambda, approx_sigma_lateral(lambda, f, unfocused_sigma, n))

#sigma of a focused beam as a function of displacement from the beam waist
displaced_sigma(lambda, focused_sigma, displacement) = sqrt(focused_sigma^2*(1+(lambda*displacement/(pi*focused_sigma^2))^2))

#This doesn't match the calculation on the OpenSPIM website.  Seems theirs is incorrect.
fwhm_stats(std) = 2*sqrt(2*log(2)) * std
fwhm_optics(sigma) = fwhm_stats(sigma/2)

#Currently aren't using these
abbe_lat(lambda, NA) = lambda/(2*NA)
abbe_ax(lambda, NA) = 2*lambda/(NA^2)
