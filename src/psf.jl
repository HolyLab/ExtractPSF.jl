#An N-dimensional PSF is modeled as a product of N gaussians oriented orthogonally
#The magnitude of the PSF at a given location in space can be found by
#indexing the Psf.
#The gaussians are centered at the origin; indexing at [0,0...] returns the maximum value
#See also: PsfFit
struct Psf{T,N} <: AbstractArray{T,N}
    sigmas::NTuple{N,T} #beam halfwidth, i.e. sigma in y = e^-(2x^2/sigma^2)
end
Psf(args...) = Psf((args...,))

function Psf(;NA=0.5, M=10.0, lambda_nm=500, um_per_pix=6.5, um_per_slice=5.0, n=1.333)
    lateral_optical_um = approx_sigma_lateral(lambda_nm/1000, NA)
    axial_optical_um = approx_sigma_axial(lambda_nm/1000, NA, n)
    lateral_pixel = lateral_optical_um/(um_per_pix/M)
    axial_pixel = axial_optical_um/um_per_slice
    return Psf((lateral_pixel, lateral_pixel, axial_pixel))
end

size_optics(psf::Psf{T,N}) where {T,N} = psf.sigmas
size(psf::Psf{T,N}) where {T,N} = size_optics(psf)
size_stats(psf::Psf{T,N}) where {T,N} = size_optics(psf) ./ 2

axes(psf::Psf{T,N}) where {T,N} = ntuple(i->nothing, N)
Base.show(io::IO, psf::Psf{T,N}) where {T,N} = print(io, "$N-dimensional Psf with size $(size(psf)) (halfwidths of Gaussian approximation)\n")
Base.show(io::IO, ::MIME"text/plain", psf::Psf) = show(io, psf)

getindex(psf::Psf{T,N}, args...) where {T,N} = getindex(psf, T.(args)...)
getindex(psf::Psf{T,N}, I::Vararg{T,N}) where {T,N} = prod(eval_gaussian.((I...,), size(psf) ./ 2))
setindex!(psf::Psf, args...) = error("setindex! is not appropriate for a Psf")

fwhm(psf::Psf, dim::Int) = fwhm_optics(size(psf)[dim])
fwhm_stats(psf::Psf, dim::Int) = fwhm_stats(size(psf)[dim])

lightsheet_psf(det_psf::Psf{T,3}, NA_illum, lambda_illum_nm, um_per_slice =1) where {T} =
    lightsheet_psf(det_psf, approx_sigma_lateral(lambda_illum_nm/1000, NA_illum)/um_per_slice)

function lightsheet_psf(det_psf::Psf{T,3}, sheet_sigma) where {T}
    dax = size(det_psf)[3]
    combinedax = lightsheet_axial_combined_psf(dax, sheet_sigma)
    return Psf((size(det_psf)[1:2]..., combinedax))
end

#PsfFit contains the result of fitting a point spread function to an image img 
#See also: Psf
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

#Fits a gaussian along each dimension independently
function PsfFit(img::OffsetArray{Timg,N}, bias;
                    guess_sigmas=ones(N),
                    guess_scales=ones(N),
                    precision=Float64) where {Timg,N}
    stds = zeros(precision,N)
    shifts, scales, quality = copy(stds), copy(stds), copy(stds)
    pimg = parent(img)
    for i = 1:N
        reduced = sum_otherdims(pimg, i, precision)
        newbias = precision(bias)*(div(length(img), size(img,i))) #sum the pixelwise bias in this slice
        curdata = reshape(reduced, length(reduced))
        curdata = OffsetArray(curdata, axes(img)[i])
        stds[i], shifts[i], scales[i], fval = fit_gaussian_1d(curdata, newbias;
                                                             guess_std=guess_sigmas[i] / 2,
                                                             guess_scale=guess_scales[i])
        quality[i] = 1.0-fval/sum(curdata.^2)
    end
    sigmas = stds .* 2 #convert from std to optics halfwidth convention
    psf = Psf((sigmas...,))
    fit = PsfFit(img, eltype(img)(bias), psf, (shifts...,), (scales...,), (quality...,))
    return fit
end
