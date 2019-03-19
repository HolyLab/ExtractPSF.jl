#generate synthetic psfs

function synth_psf_img(psf::Psf{T,N}, imgradii::NTuple{N,Int}; shift=(zeros(N)...,), scale=one(T), bias=zero(T)) where {T,N}
    imginds = roiranges(imgradii)
    #imgdat = zeros(T, 2 .* imgradii .+ 1)
    imgdat = fill(T(bias), (2 .* imgradii .+ 1))
    img = OffsetArray(imgdat, imginds)
    synth_psf_img!(img, psf; shift=shift, scale=scale)
end

mayberound(t::Val{T}, v) where {T<:Integer} = round(T, v)
mayberound(t, v) = v

#adds a psf onto existing input img
function synth_psf_img!(img::OffsetArray{TA,N}, psf::Psf{T,N}; shift=(zeros(N)...,), scale=one(T)) where {TA,T,N}
    sigmas = size(psf)
    for I in CartesianIndices(img)
        img[I] += mayberound(Val(TA), scale*psf[(Tuple(I).+shift)...])
    end
    return img
end

rng_subtr(r::UnitRange, v::Int) = (first(r)-v):(last(r)-v)
rng_subtr(r, v::Int) = r-v

#because subtracting a tuple of ints from a tuple of ranges doesn't seem to work...
tuple_subtr(t1::Tuple, t2::Tuple) = (rng_subtr(first(t1), first(t2)), tuple_subtr(Base.tail(t1), Base.tail(t2))...)
tuple_subtr(t1::Tuple{}, t2::Tuple{}) = ()

#noise will be uniformly distributed, +- noise_max
function fake_bead_sample(img_size::NTuple{3,Int}, bead_sigma, bead_shift, bead_ctrs::AbstractVector, bead_max, noise_max, bias)
    roi_radii = 2 .* ceil.(Int, bead_sigma)
    psf0 = Psf(bead_sigma)
    fullimg = zeros(UInt16, img_size...)
    for i in eachindex(fullimg)
        fullimg[i] = round(eltype(fullimg), bias + (rand()-0.5)*2*noise_max)
    end
    rngs0 = ExtractPSF.roiranges(roi_radii) #indices at which to index psf if its roi is contained in image
    for c in bead_ctrs
        crngs = ExtractPSF.roiranges(c, roi_radii) #parent indices we'd like to use
        valid_rngs = intersect.(crngs, axes(fullimg)) #valid indices in parent array
        vshifted = tuple_subtr(valid_rngs, c)
        psf_rngs = intersect.(rngs0, vshifted) #evaluation indices for psf
        imgv = view(fullimg, valid_rngs...)
        oa = OffsetArray(imgv, psf_rngs)
        ExtractPSF.synth_psf_img!(oa, psf0; shift=bead_shift, scale=bead_max)
    end
    return fullimg
end
