using ExtractPSF
using Test

_sigma = (5.0,10.0,15.0)
bias = 100
_scale = 2.0
_shift = (1.4,0,0)
psf0 = Psf(_sigma)
img0 = ExtractPSF.synth_psf_img(psf0, (30,30,30); shift=_shift, scale=_scale, bias=bias)
#add some noise
noise = OffsetArray((rand(size(img0)...).-0.5).*0.5, axes(img0))
img = img0.+noise

#fit psf
fit = PsfFit(img, bias)
psff = psf(fit)
shiftf = shift(fit)
q = quality(fit)

@test all(q.>0.9999)
for (w0,w1) in zip(size(psf0), size(psff))
    @test abs(w0-w1) < 0.3
end
for (s0,s1) in zip(_shift, shift(fit))
    @test abs(s0-s1) < 0.3
end
