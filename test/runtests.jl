using QWNNM
using Test

nSig  = 10
Img_O = loadImg("lena.png")
Img_N = loadImg("lenaN.png")

@test psnr(addNoise(Img_O, 10), Img_O) ≈ psnr(Img_N, Img_O)


