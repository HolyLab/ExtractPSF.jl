using ExtractPSF, OffsetArrays
using Test


bead_sigma = (1.0,2.0,3.0)
roi_radii = 2 .* ceil.(Int, bead_sigma)
bias = 100
bead_shift = (0,0,0)
noise_max = 3
bead_max = 100.0 #bead peak intensity
imgsz = (50,50,20)

#Case: first bead's roi isn't fully contained, so it should be excluded
bead_ctrs = [(9,10,6), (30, 10, 7)]

fullimg = ExtractPSF.fake_bead_sample(imgsz, bead_sigma, bead_shift, bead_ctrs, bead_max, noise_max, bias)

#rois is an array of OffsetArrays, ctrs is an array of coords (x,y,z)
rois, ctrs = find_beads(fullimg, roi_radii, bias)

@test length(rois) == 1
@test ctrs[1] == bead_ctrs[2]
@test all(rois[1] .== OffsetArray(view(fullimg, 28:32, 6:14, 1:13), -2:2, -4:4, -6:6))

#Case: bead's peak is at saturation, so should be excluded
fullimg[bead_ctrs[2]...] = typemax(eltype(fullimg))
rois, ctrs = find_beads(fullimg, roi_radii, bias)
@test length(rois) == 0

#Case: bead ROIs overlap, so both should be excluded
bead_ctrs = [(26,10,6), (30, 10, 7)]
fullimg = ExtractPSF.fake_bead_sample(imgsz, bead_sigma, bead_shift, bead_ctrs, bead_max, noise_max, bias)
rois, ctrs = find_beads(fullimg, roi_radii, bias)
@test length(rois) == 0
