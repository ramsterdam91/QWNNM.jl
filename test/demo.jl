using QWNNM
using ImageView

function loadImg(path)
    #from [0.,1.] to [0,256] 
    Img = 255 * float32.(channelview(load(path)))
    #Switch 1,3 Axis
    Img = permutedims(Img,(2, 3, 1))
    #Extract RGB only
    Img = Img[:,:,1:3]
end

function beImg(Img::Array{Float32,3})
    Img = permutedims(Img,(3, 1, 2))./255
    Img = colorview(RGB, Img)
end

nSig  = 10
Img_O = loadImg("lena.png")
Img_N = loadImg("lenaN.png")

Par = QWNNM.autopar(nSig)

Img_E = QWNNM.denoising(Img_O, Img_N, Par)

imshow(beImg(Img_E))