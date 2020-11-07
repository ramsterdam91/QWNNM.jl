function WNNM(Y::Array{Float32,3}, C::Float32, Sig::Array{Float32,1}, m::Array{Float32,3})

    F = svd(F2cplx(Y))

    lenSig = length(F.S)

    patNum = size(Y)[2]

    mean_Sig_arr = √(sum(Sig[:].^2) / length(Sig[:]))

    TempC  = Float32(C * √(patNum) * 2 * mean_Sig_arr^2)

    X_arr = ClosedQWNNM(F.S, TempC,eps())

    X = F.U * Diagonal(X_arr) * F.Vt

    return cplx2F(X) .+ m
end

function ClosedQWNNM(SigmaY::Array{Float32,1}, c::Float32, our_eps)

    temp = (SigmaY .- our_eps) .^ 2 - 4 * ((-1) * our_eps .* SigmaY .+ c)
    ind = Int32[]
    SigmaX = Float32[]
    for i in 1:length(temp)
        if temp[i] > 0
            push!(ind, i)
            push!(SigmaX, max(SigmaY[i] - our_eps + √(temp[i]) , 0) / 2) 
        else 
            push!(SigmaX, 0)
        end
    end
    
    return SigmaX
end

function F2cplx(M::Array{Float32,3})::Array{Complex{Float32},2}
    (h, w) = size(M)

    cM = adjoint(M)

    [cM[1:2*h, 1:w]; cM[1:2*h, w+1:2*w]]
end

function adjoint(A::Array{Float32,3})::Array{Complex{Float32},2}
    W = zeros(size(A)[1:2])
    X = A[:, :, 1]
    Y = A[:, :, 2]
    Z = A[:, :, 3]
    imm = Complex{Float32}(im)
    vcat(hcat(W + imm .* X, Y + imm .* Z),hcat(-Y + imm .* Z, W - imm .* X))
end

function cplx2F(M)
    (h, w) = size(M)

    cM = [M[1:h ÷ 2,1:w] M[1+h ÷ 2:h,1 : w]]

    unadjoint(cM)
end
function unadjoint(A)
    (r, c) = size(A)
    r2 = r ÷ 2
    c2 = c ÷ 2
    A1 = A[     1 : r2,      1 : c2]
    A2 = A[     1 : r2, c2 + 1 : c ]
    A3 = A[r2 + 1 : r,       1 : c2]
    A4 = A[r2 + 1 : r,  c2 + 1 : c ]
    C = zeros(Float32, r2, c2, 3)
    C[:,:,1] = Float32.(imag(A1 - A4))
    C[:,:,2] = Float32.(real(A2 - A3))
    C[:,:,3] = Float32.(imag(A2 + A3) )
    C./2
end