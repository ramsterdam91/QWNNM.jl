#input Img_Noising and Img_Origin
#return the PSNR

function psnr(Img1::Array{Float32,3}, Img2::Array{Float32,3}) :: Float64

    if size(Img1) != size(Img2)
        println("Images should be in the same size")
        return
    end

    (N,M) = size(Img1)
  
    10 * log10(255^2 / mean((Img1 - Img2).^2))
    
end