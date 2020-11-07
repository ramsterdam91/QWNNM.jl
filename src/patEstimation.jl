#slowwwwwwwwwwwwwwwwwww
#Curent iteration Patches `CurPat`, Noise level of Patches `Sigma_arr`
function patEstimation(NL_mat::Array{Int32,2}, Self_arr::Array{Int32,2}, Sigma_arr::Array{Float32,3}, CurPat::Array{Float32,3}, Par::PAR)
    Depth = size(CurPat)[3]

    EPat   = zeros(size(CurPat))
    W      = zeros(size(CurPat))
    #mean
    M_Temp = zeros(size(CurPat))

    for i = 1 : length(Self_arr)
        Temp    =   CurPat[:, NL_mat[:, i], :]
        M_Temp  =   mean(Temp, dims = 2)
        Temp    =   Temp .- M_Temp
        E_Temp 	=   WNNM(Temp, Par.c, Sigma_arr[1, Self_arr[i], :], M_Temp)

        EPat[:,NL_mat[:, i],:]    += E_Temp
        W[:,NL_mat[:, i], :]      += ones(Par.patsize^2, size(NL_mat[:, i])[1], Depth)
    end
    
    (EPat, W)
end