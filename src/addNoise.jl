function addNoise(Img::Array{Float32,3}, nSig::Int64)

    (N, M) = size(Img)

    Img_N = Img .+ nSig * randn(N, M)

end