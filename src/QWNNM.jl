__precompile__

module QWNNM

    using Images
    using Random
    using Statistics
    using LinearAlgebra
    using Quaternions
    using GenericSVD

    include("addNoise.jl")
    include("psnr.jl")
    include("denoising.jl")
    
end
