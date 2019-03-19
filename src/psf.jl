#gaussian approximation to airy disk
#taken from paraxial approximations derived in:
#Gaussian approximations of fluorescence microscopepoint-spread function models
#Bo Zhang, Josiane Zerubia, and Jean-Christophe Olivo-Marin
#1 April 2007 Vol. 46, No. 10 APPLIED OPTICS
#gauss_lat(lambda, NA) = sqrt(2)/((2*pi)*NA/lambda)
gauss_lat(lambda, NA) = sqrt(2)/(((2*pi)/lambda)*NA)
gauss_ax(lambda, NA, n) = 2*sqrt(6)*n/(2*pi/lambda*NA^2)

#adapted from wikipedia and paper above
#na2fnum(NA, n) = n/(2*NA)
#fnum2na(fnum, n) = (n/2) * (1/fnum)
#gauss_lat2(lambda, NA, n) = 0.44/n * lambda * na2fnum(NA, n) #approx for 2D case (almost same as 3D)

#assuming a lens has a big enough aperture to accept
#a collimated input of specified diameter, what is the NA of the system?
#function lens_na(collimated_diameter, bfl, f, n)
#    r = collimated_diameter/2
#    return n * r/sqrt(r^2 + f^2)
#end

abbe_lat(lambda, NA) = lambda/(2*NA)
abbe_ax(lambda, NA) = 2*lambda/(NA^2)

gaussian_fwhm(sigma) = 2*sqrt(2*log(2))*sigma
#1/e^2 beamwidth (full) from fwhm for a gaussian beam
fwhm2beamwidth(fwhm) = sqrt(2)*fwhm/sqrt(log(2))
gaussian2beamwidth(sigma) = fwhm2beamwidth(gaussian_fwhm(sigma))

function rayleigh_length(NA, lambda)
    bw = gaussian2beamwidth(gauss_lat(lambda, NA))
    return pi*bw^2/lambda
end


struct Psf{T,N} <: AbstractArray{T,N}
    sigmas::NTuple{N,T}
end
Psf(args...) = Psf((args...,))

function Psf(;NA=0.5, M=10.0, lambda_nm=500, um_per_pix=6.5, um_per_slice=5.0, n=1.333)
    lateral_optical_um = gauss_lat(lambda_nm/1000, NA)
    axial_optical_um = gauss_ax(lambda_nm/1000, NA, n)
    lateral_pixel = lateral_optical_um/(um_per_pix/M)
    axial_pixel = axial_optical_um/um_per_slice
    return Psf((lateral_pixel, lateral_pixel, axial_pixel))
end

size(psf::Psf{T,N}) where {T,N} = psf.sigmas

axes(psf::Psf{T,N}) where {T,N} = ntuple(i->nothing, N)
Base.show(io::IO, psf::Psf{T,N}) where {T,N} = print(io, "$N-dimensional Psf with size $(size(psf)) (stds of Gaussian approximation)\n")
Base.show(io::IO, ::MIME"text/plain", psf::Psf) = show(io, psf)


getindex(psf::Psf{T,N}, args...) where {T,N} = getindex(psf, T.(args)...)
getindex(psf::Psf{T,N}, I::Vararg{T,N}) where {T,N} = prod(eval_gaussian.((I...,), size(psf)))
setindex!(psf::Psf, args...) = error("setindex! is not appropriate for a Psf")

gaussian_lightsheet_psf(det_psf::Psf{T,3}, NA_illum, lambda_illum_nm, um_per_slice =1) where {T} =
    #gaussian_lightsheet_psf(det_psf, abbe_lat(lambda_illum_nm/1000, NA_illum)/um_per_slice)
    gaussian_lightsheet_psf(det_psf, gauss_lat(lambda_illum_nm/1000, NA_illum)/um_per_slice)
function gaussian_lightsheet_psf(det_psf::Psf{T,3}, sheet_sigma) where {T}
    dax = size(det_psf)[3]
    combinedax = sqrt((dax^2*sheet_sigma^2)/(dax^2 + sheet_sigma^2))
    return Psf((size(det_psf)[1:2]..., combinedax))
end

#1D gaussian centered at zero, peak height of 1
eval_gaussian(x, sigma) = exp(-((x^2)/(2*sigma^2)))

struct PsfFit{Tfit,N,Timg}
    img::OffsetArray{Timg,N} #image array to which we fit the psf
    bias::Timg #black level of image (100 for PCO cameras)
    psf::Psf{Tfit,N}
    shift::NTuple{N,Tfit} #amount by which fitted psf is shifted relative to the center of the img array
    scale::NTuple{N,Tfit} #amount by which fitted psf is scaled (scale of 1 means peak height is 1)
    quality::NTuple{N,Tfit} #fractions of variance explained by the gaussian fits in each dimension
end

psf(fit::PsfFit) = fit.psf
quality(fit::PsfFit) = fit.quality
size(fit::PsfFit) = size(psf(fit))
shift(fit::PsfFit) = fit.shift
scale(fit::PsfFit) = fit.scale
image(fit::PsfFit) = fit.img

function gaussian_obj(sigma, shift, scale, bias, data::OffsetArray{T,1}) where {T}
    if sigma > 0 && scale > 0
        p = zero(sigma)
        for i in first(axes(data))
            #squared error of fitted gaussian and data, evaluated at each index of the data array
            p += (data[i]-(bias+scale*eval_gaussian(i+shift, sigma)))^2
        end
        return p
    else #infeasible sigma, scale
        m = min(sigma, scale)
        return sum(data.^2)*(1+m^2)
    end
end

gaussdecode(m) = (m[1], m[2], m[3])

#sigma, shift, and scale packed into params array
gaussian_obj(params, bias, data) = gaussian_obj(gaussdecode(params)..., bias, data)

function sum_otherdims(img, dim::Int, precision)
    sz = ntuple(d->d!=dim ? 1 : size(img,dim), Val(ndims(img)))
    reduced = Array{precision}(undef, sz...)
    sum_otherdims!(reduced, img)
end

function sum_otherdims!(reduced, img)
    fill!(reduced, 0)
    rmax = last(CartesianIndices(reduced))
    for I in CartesianIndices(img)
        reduced[min(rmax, I)] += img[I]
    end
    return reduced
end

function PsfFit(img::OffsetArray{Timg,N}, bias;
                    guess_sigmas=ones(N),
                    guess_scales=ones(N),
                    precision=Float64) where {Timg,N}
    sigmas = zeros(precision,N)
    shifts, scales, quality = copy(sigmas), copy(sigmas), copy(sigmas)
    pimg = parent(img)
    for i = 1:N
        reduced = sum_otherdims(pimg, i, precision)
        newbias = bias*(div(length(img), size(img,i))) #sum the pixelwise bias in this slice
        curdata = reshape(reduced, length(reduced))
        curdata = OffsetArray(curdata, axes(img)[i])
        gobj = x->gaussian_obj(x, newbias, curdata)
        firstguess = precision.([guess_sigmas[i]; 0; guess_scales[i]])
        result = optimize(gobj, firstguess, LBFGS(); autodiff=:forward)
        m = Optim.minimizer(result)
        fval = Optim.minimum(result)
        quality[i] = 1.0-fval/sum(curdata.^2)
        sigmas[i], shifts[i], scales[i] = gaussdecode(m)
    end
    psf = Psf((sigmas...,))
    fit = PsfFit(img, eltype(img)(bias), psf, (shifts...,), (scales...,), (quality...,))
    return fit
end
