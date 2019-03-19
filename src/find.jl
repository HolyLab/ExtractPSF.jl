##check whether two ranges intersect
#function intersectrng(r1::UnitRange, r2::UnitRange)
#    i = intersect(r1)
#    return first(i)<=last(i)
#end
##check whether all corresponding ranges in two lists of ranges intersect
#intersectrng(t1::Tuple, t2::Tuple) = intersectrng(first(t1), first(t2)) && intersectrng(Base.tail(t1), Base.tail(t2))
#intersectrng(t1::Tuple{}, t2::Tuple{}) = true

crop_k_inds(img::AbstractArray, k::Int) = crop_k_inds(axes(img), k)
crop_k_inds(inds::Tuple, k::Int) = ((first(first(inds))+k):(last(first(inds))-k), crop_k_inds(Base.tail(inds), k)...)
crop_k_inds(inds::Tuple{}, k::Int) = ()

#negranges(t::Tuple) = (-first(t):first(t), negranges(Base.tail(t))...)
roiranges(r::Tuple) = roiranges(ntuple(i->zero(first(r)),length(r)), r) #assume centered at origin
roiranges(ctr::Tuple, r::Tuple) = ((first(ctr)-first(r)):(first(ctr)+first(r)), roiranges(Base.tail(ctr), Base.tail(r))...)
roiranges(ctr::Tuple{}, r::Tuple{}) = ()

#return iterator over all indices adjacent to I
neighbor_inds(I::CartesianIndex) = filter(x->x!=I, CartesianIndices(roiranges(I.I, ntuple(i->1, length(I)))))

function find_local_maxima(img::AbstractArray, thresh)
    inds = CartesianIndices(crop_k_inds(img, 1))
    maxima = typeof(first(inds).I)[]
    for I in inds
        v = img[I]
        if v >= thresh
            is_maximum = true
            for i in neighbor_inds(I)
                if v <= img[i]
                    is_maximum = false
                    break
                end
            end
            if is_maximum
                push!(maxima, I.I)
            end
        end
    end
    return maxima
end

#for each local intensity maximum assign an integer index and paint a canvas with it
#if while painting the canvas we encounter an already-painted pixel then this roi overlaps
#...so we exclude it from the list of found beads
function exclude_overlapping(roi_ctrs, roi_radii, imgsz)
    canvas = zeros(Int, imgsz)
    keep = trues(length(roi_ctrs))
    for (i,c) in enumerate(roi_ctrs)
        #printed = false
        for I in CartesianIndices(roiranges(c, roi_radii))
            if checkbounds(Bool, canvas, I) #could be more efficient
                v = canvas[I]
                if v!=0
                    #if !printed
                        #print("overlap between roi centered at $c and roi centered at $(roi_ctrs[v])\n")
                        #print("radius was $roi_radii \n")
                        #printed=true
                    #end
                    keep[i] = false
                else
                    canvas[I] = i
                end
            else #the roi extends outside of the image
                keep[i] = false
            end
        end
    end
    return keep
end

#find beads that do not overlap within radiii roi_radii and return a list of rois (views into parent image)
#along with a list of coordinate centers in parent image
function find_beads(fullimg::AbstractArray{T,N}, roi_radii::NTuple{N,Int}, bias; noise_stds=5, sat_val=typemax(T)) where {T,N}
    #estimate noise level in black pixels, assuming at least 80% of pixels are "black"
    psorted = sort(reshape(fullimg, length(fullimg)))
    noise = std(view(psorted, 1:round(Int, 0.8*length(psorted))))
    #blur the image to eliminate multiple local maxima from nearby pixels of same bead
    kern = KernelFactors.gaussian(div.(roi_radii,2))
    kernmax = maximum(first(kern))
    imgf = imfilter(fullimg, kern)
    #find all local intensity maxima that are above a noise threshold (set by estimated noise)
    thresh = bias+kernmax*noise*noise_stds
    maxima = find_local_maxima(imgf, thresh)
    keep = exclude_overlapping(maxima, roi_radii, size(fullimg))
    if all(keep.==false)
        @warn("No valid bead ROIs found.  Common reasons for this are that the bead sample was too dense, or that the image is too noisy.")
        return [],[]
    else
        mkeep = maxima[keep]
        firstroi = view(fullimg, roiranges(mkeep[1], roi_radii)...)
        firstroi = OffsetArray(firstroi, roiranges(roi_radii)...)
        rois = typeof(firstroi)[]
        keep2 = trues(length(mkeep))
        for (i,c) in enumerate(mkeep)
            roi = view(fullimg, roiranges(c, roi_radii)...)
            roi = OffsetArray(roi, roiranges(roi_radii)...)
            if all(roi.<sat_val)
                push!(rois, roi)
            else
                #print("Excluding ROI centered at $c due to saturation\n")
                keep2[i] = false
            end
        end
        if all(keep2.==false)
            @warn("After removing ROIs with saturated pixels no valid bead ROIs remained.  Use less intense illumination or a shorter exposure time.")
            return [], []
        end
        return rois, mkeep[keep2]
    end
end
